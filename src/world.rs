use std::ops::Mul;
use crate::BVH::BVH;

use crate::camera::{BufferItem, Camera};
use crate::geometric::{Point, Ray, Triangle, Vector};
use crate::light::{Color, Light};
use crate::object::Object;

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

    fn coloring(&self, ray: &Ray, triangle: Option<((&Triangle, &Object, usize), f64, Point)>) -> Color {
        if let Some(((triangle, object, index), _, point)) = triangle {
            let mut color = Color::black();
            for light in &self.lights {
                let light_ray = Ray {
                    start: light.position.clone(),
                    direction: (Vector::from(&point) - Vector::from(&light.position)).normalize(),
                };
                let triangle_index = self.bvh.as_ref().unwrap().intersect(&light_ray);
                let shadowed = if let Some(triangle_index) = triangle_index {
                    triangle_index != index
                } else {
                    false
                };
                let normal = get_normal(&point, &triangle, NormMixMode::VertexDistanceReverseFaceAverage);
                color = color + light.phong(&point, normal, ray.direction, &object.properties, shadowed);
            }
            color
        } else {
            Color::default()
        }
    }

    pub fn render(&mut self) {
        // count faces
        let mut count = 0;
        for object in &self.objects {
            count += object.triangles.len();
        }
        println!("Rendering {} objects and {} faces", self.objects.len(), count);
        
        // calculate normals
        for object in &mut self.objects {
            object.calc_triangle_norms();
            object.calc_point_norms();
        }

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
        for object in &mut self.objects {
            object.calc_triangle_norms();
        }
        
        // project to camera buffer
        let mut i = 1;
        let mut triangle_object_s = vec![];
        for object in &self.objects {
            for triangle in &object.triangles {
                self.camera.project(triangle, i);
                triangle_object_s.push((triangle, object, i - 1));
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
        for y in 0..self.camera.picture.height {
            for x in 0..self.camera.picture.width {
                let ray = self.camera.get_ray(x, y);
                let BufferItem {index, depth, point} = self.camera.buffer[(x as usize, y as usize)].clone();
                let color = self.coloring(&ray, if index == 0 { None } else { Some((triangle_object_s[index - 1], depth, point)) });
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