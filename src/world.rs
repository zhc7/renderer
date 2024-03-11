use std::collections::BinaryHeap;
use std::ops::Mul;
use crate::BVH::BVH;

use crate::camera::{BufferItem, Camera};
use crate::geometric::{Point, Ray, Triangle, Vector};
use crate::light::{Color, Light};
use crate::object::Object;


use tqdm::tqdm;

#[derive(Clone)]
pub struct Matrix {
    pub m: [[f64; 4]; 4],
}

impl Matrix {
    fn shift(vector: &Vector) -> Matrix {
        Matrix {
            m: [
                [1.0, 0.0, 0.0, vector.x],
                [0.0, 1.0, 0.0, vector.y],
                [0.0, 0.0, 1.0, vector.z],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Matrix) -> Matrix {
        let mut result = Matrix {
            m: [[0.0; 4]; 4],
        };
        for i in 0..4 {
            for j in 0..4 {
                result.m[i][j] = self.m[i][0] * rhs.m[0][j] +
                    self.m[i][1] * rhs.m[1][j] +
                    self.m[i][2] * rhs.m[2][j] +
                    self.m[i][3] * rhs.m[3][j];
            }
        }
        result
    }
}

pub struct World {
    pub objects: Vec<Object>,
    pub lights: Vec<Light>,
    pub camera: Camera,
    bvh: Option<BVH>,
}

pub enum NormMixMode {
    Flat,
    VertexDistanceReverseAverage,
    VertexDistanceReverseFaceAverage,
    Phong,
}

fn get_normal(point: &Point, triangle: &Triangle, mode: NormMixMode) -> Vector {
    match mode {
        NormMixMode::Flat => triangle.normal.unwrap(),
        NormMixMode::VertexDistanceReverseAverage | NormMixMode::VertexDistanceReverseFaceAverage => {
            // average normal of three vertices
            let mut normal = Vector::zero();
            let mut weights = 0.0;
            for i in 0..3 {
                let ver = &triangle[i];
                let w = 1. / point.distance(&ver);
                normal = normal + ver.normal().unwrap() * w;
                weights += w;
            }
            normal = normal / weights;
            match mode {
                NormMixMode::VertexDistanceReverseAverage => normal,
                NormMixMode::VertexDistanceReverseFaceAverage => (normal + triangle.normal.unwrap()) / 2.0,
                _ => panic!(),
            }
        },
        NormMixMode::Phong => {
            // phong interpolation
            let mut normal = Vector::zero();
            let mut weights = 0.0;
            for i in 0..3 {
                let a = &triangle[i];
                let b = &triangle[(i + 1) % 3];
                let c = &triangle[(i + 2) % 3];
                let w = (Vector::from(b) - Vector::from(point)).cross(Vector::from(c) - Vector::from(point)).magnitude();
                normal = normal + a.normal().unwrap() * w;
                weights += w;
            }
            normal = normal / weights;
            normal
        }
    }.normalize()
}

fn refract(ray: &Ray, normal: Vector, n1: f64, n2: f64) -> (Option<Vector>, f64) {
    let cos_theta_i = -ray.direction.dot(normal);
    let sin_theta_i = (1.0 - cos_theta_i.powi(2)).sqrt();
    let sin_theta_t = n1 / n2 * sin_theta_i;
    if sin_theta_t > 1.0 {
        return (None, 1.0);
    }
    let cos_theta_t = (1.0 - sin_theta_t.powi(2)).sqrt();
    let r_parallel = (n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t);
    let r_perpendicular = (n2 * cos_theta_i - n1 * cos_theta_t) / (n2 * cos_theta_i + n1 * cos_theta_t);
    let reflect = (r_parallel.powi(2) + r_perpendicular.powi(2)) / 2.0;
    let refract_direction = (ray.direction + normal * cos_theta_i) * (n1 / n2) - normal * cos_theta_t;
    (Some(refract_direction), reflect)
}

#[derive(Clone, Copy)]
struct Status {
    // inside which object, at what depth
    index: usize,
    depth: f64,
    refractive_index: f64,
    transparency: f64,
}

impl Status {
    fn air(depth: f64) -> Status {
        Status {
            index: usize::MAX,
            depth,
            refractive_index: 1.0,
            transparency: 1.0,
        }
    }
}

impl World {
    pub fn new() -> World {
        World {
            objects: Vec::new(),
            lights: Vec::new(),
            camera: Camera::new(),
            bvh: None,
        }
    }

    pub fn add_object(&mut self, object: Object, center: Point) {
        for point in &object.points {
            point.transform(&Matrix::shift(&Vector::from(&center)));
        }
        self.objects.push(object);
    }

    pub fn add_light(&mut self, light: Light, center: Point) {
        light.position.transform(&Matrix::shift(&Vector::from(&center)));
        self.lights.push(light);
    }

    fn transform(&mut self, matrix: Matrix) {
        for object in &mut self.objects {
            for point in &mut object.points {
                point.transform(&matrix);
            }
        }
        for light in &mut self.lights {
            light.position.transform(&matrix);
        }
        self.camera.position.transform(&matrix);
    }

    fn coloring(&self, ray: &Ray, buffered: &mut BinaryHeap<BufferItem>, seq: &Vec<(&Triangle, &Object, usize)>, status: Status, mut ttl: u32) -> Color {
        if ttl == 0 {
            return Color::black();
        }
        ttl -= 1;
        if buffered.peek().is_none() {
            return Color::default();
        }
        let item = buffered.pop().unwrap();
        let (triangle, object, index) = seq[item.index];
        let point = &item.point;
        let mut color = Color::black();
        let absorbed = status.transparency.powf(item.depth);


        let normal = if triangle.size() < 100.0 {
            get_normal(&point, &triangle, NormMixMode::Phong)
        } else {
            triangle.normal.unwrap()
        };
        assert!(1.0 - normal.magnitude() < 1e-6);

        // direct lights
        if status.index == usize::MAX {
            for light in &self.lights {
                let light_ray = Ray {
                    start: light.position.clone(),
                    direction: (Vector::from(point) - Vector::from(&light.position)).normalize(),
                };
                let triangle_indices = self.bvh.as_ref().unwrap().intersect(&light_ray);

                // calculate remaining light, ignore refraction
                let mut current = 1.0;
                // we assume light is not inside anything
                let mut status = Status::air(0.0);

                for BufferItem { index: i, depth: t, point: _ } in triangle_indices.into_sorted_vec().into_iter().rev() {
                    if i == index {
                        // reached
                        break;
                    }
                    let object = seq[i].1;
                    if status.index < usize::MAX {
                        current *= object.properties.transparent.powf(t - status.depth);
                        if current == 0.0 {
                            break;
                        }
                        status = Status::air(t);
                    } else {
                        status = Status {
                            index: i,
                            depth: t,
                            refractive_index: object.properties.refractive_index,
                            transparency: object.properties.transparent,
                        };
                    }
                }

                color = color + light.phong(&point, normal, ray.direction, &object.properties, current);
            }
        }


        let normal = if status.index < usize::MAX {
            // we are at the inner side of the object
            -normal
        } else {
            normal
        };
        // Although the triangle is facing the camera, this point might be behind the camera because of interpolation,
        // we don't calculate its reflection and transparency as it may lead to bigger mistakes.
        let cos_theta_i = -ray.direction.dot(normal);
        if ttl > 0 && cos_theta_i > 0.0 {
            // reflected color
            let tri_angle = ray.direction.dot(triangle.normal.unwrap());
            assert_eq!(tri_angle < 0.0, status.index == usize::MAX);
            let reflection_ray = Ray {
                start: point.clone(),
                direction: ray.direction + normal * 2.0 * cos_theta_i,
            };
            let mut reflect_hit = self.bvh.as_ref().unwrap().intersect(&reflection_ray);
            if reflect_hit.len() % 2 == (status.index == usize::MAX) as usize {
                // those do not meet this requirement are also because of the interpolation
                // when in air, should hit even times
                return color * absorbed;
            }
            // use fresnel equation to calculate the reflectivity
            let n1: f64 = status.refractive_index;
            let n2: f64 = if status.index == usize::MAX {
                object.properties.refractive_index
            } else {
                1.0
            };
            let (refract_direction, reflect) = refract(ray, normal, n1, n2);
            let refract = 1.0 - reflect;
            color += self.coloring(&reflection_ray, &mut reflect_hit, seq, Status {
                index: if status.index == usize::MAX {
                    // from air, to air
                    usize::MAX
                } else {
                    index
                },
                depth: item.depth + status.depth,
                // still inside the same medium
                refractive_index: status.refractive_index,
                transparency: status.transparency,
            }, ttl) * reflect;

            // transparent color
            if object.properties.transparent == 0.0 || refract_direction.is_none() {
                return color * absorbed;
            }
            let refract_ray = Ray {
                start: point.clone(),
                direction: refract_direction.unwrap(),
            };
            let mut refract_hit = self.bvh.as_ref().unwrap().intersect(&refract_ray);
            if refract_hit.len() % 2 == (status.index < usize::MAX) as usize {
                // those do not meet this requirement are also because of the interpolation
                // when in air, should hit odd times
                return color * absorbed;
            }
            color += self.coloring(&refract_ray, &mut refract_hit, seq, Status {
                index: if status.index == usize::MAX {
                    // from air, to object
                    index
                } else {
                    usize::MAX
                },
                depth: item.depth + status.depth,
                refractive_index: n2,
                transparency: object.properties.transparent,
            }, ttl) * refract;
        }

        color * absorbed
    }

    pub fn render(&mut self) {
        // count faces
        let mut count = 0;
        for object in &self.objects {
            count += object.triangles.len();
        }
        println!("Rendering {} objects and {} faces", self.objects.len(), count);

        // transform to camera view
        let mat_shift = Matrix::shift(&-Vector::from(&self.camera.position));
        let mat_rotate = Matrix {
            m: [
                [self.camera.right.x, self.camera.right.y, self.camera.right.z, 0.0],
                [self.camera.up.x, self.camera.up.y, self.camera.up.z, 0.0],
                [self.camera.direction.x, self.camera.direction.y, self.camera.direction.z, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };
        let mat_camera = mat_rotate * mat_shift;
        self.transform(mat_camera);
        self.camera.direction = Vector::new(0.0, 0.0, 1.0);
        self.camera.up = Vector::new(0.0, 1.0, 0.0);
        self.camera.right = Vector::new(1.0, 0.0, 0.0);

        // calculate normals
        for object in &mut self.objects {
            object.calc_triangle_norms();
            object.calc_point_norms();
        }

        // project to camera buffer
        let mut i = 0;
        let mut triangle_object_s = vec![];
        for object in &self.objects {
            for triangle in &object.triangles {
                self.camera.project(triangle, i);
                triangle_object_s.push((triangle, object, i));
                i += 1;
            }
        }

        // create bvh
        let mut triangles = vec![];
        for object in &self.objects {
            for triangle in &object.triangles {
                triangles.push(triangle);
            }
        }
        let mut bvh = BVH::new();
        bvh.build(&triangles);
        self.bvh = Some(bvh);

        // render
        for y in tqdm(0..self.camera.picture.height) {
            for x in 0..self.camera.picture.width {
                let ray = self.camera.get_ray(x, y);
                let mut buffered = self.camera.buffer[(x as usize, y as usize)].clone();
                let color = self.coloring(&ray, &mut buffered, &triangle_object_s, Status::air(0.0), 10);
                self.camera.picture[(x as usize, y as usize)] = color;
                // self.camera.picture[(x as usize, y as usize)] = if self.camera.buffer[(x as usize, y as usize)].0 != usize::default() {
                //     Color::black()
                // } else {
                //     Color::default()
                // };
            }
        }
    }
}