use std::ops::Mul;

use crate::camera::Camera;
use crate::geometric::{Ray, Vector};
use crate::light::{Color, Light};
use crate::object::Object;

#[derive(Clone)]
pub struct Matrix {
    pub m: [[f64; 4]; 4],
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

pub struct World<'a> {
    pub objects: Vec<Object<'a>>,
    pub lights: Vec<Light>,
    pub camera: Camera,
}

impl<'a> World<'a> {
    pub fn new() -> World<'a> {
        World {
            objects: Vec::new(),
            lights: Vec::new(),
            camera: Camera::new(),
        }
    }

    fn add_object(&mut self, object: Object<'a>) {
        self.objects.push(object);
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
        self.camera.direction.transform(&matrix);
        self.camera.up.transform(&matrix);
        self.camera.right.transform(&matrix);
    }

    fn trace(&self, ray: &Ray) -> Color {
        let mut nearest = None;
        let mut nearest_distance = f64::INFINITY;
        for object in &self.objects {
            for triangle in &object.triangles {
                if let Some((distance, point)) = triangle.intersect(ray) {
                    if distance < nearest_distance {
                        nearest_distance = distance;
                        nearest = Some((point, object, triangle));
                    }
                }
            }
        }
        if let Some((point, object, triangle)) = nearest {
            let mut color = Color::default();
            // average normal of three vertices
            let mut normal = Vector::new();
            let mut weights = 0.0;
            for i in 0..3 {
                let ver = &triangle[i];
                let distance = point.distance(ver);
                normal = normal + point.normal.unwrap() * distance;
                weights += distance;
            }
            normal = normal / weights;
            for light in &self.lights {
                color = color + light.phong(&point, normal, ray.direction, &object.properties);
            }
            color
        } else {
            Color::default()
        }
    }

    pub fn render(&mut self) {
        // transform to camera view
        let mat_shift = Matrix {
            m: [
                [1.0, 0.0, 0.0, -self.camera.position.x],
                [0.0, 1.0, 0.0, -self.camera.position.y],
                [0.0, 0.0, 1.0, -self.camera.position.z],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };
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

        // render
        for y in 0..self.camera.picture.height {
            for x in 0..self.camera.picture.width {
                let ray = self.camera.get_ray(x, y);
                let color = self.trace(&ray);
                self.camera.picture[(x as usize, y as usize)] = color;
            }
        }
    }
}