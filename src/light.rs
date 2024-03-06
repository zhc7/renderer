use std::ops::Mul;
use std::ops::Add;
use crate::geometric::{Point, Vector};
use crate::object::Properties;


#[derive(Clone, Copy)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

pub struct Light {
    pub color: Color,
    pub position: Point,
}

impl Default for Color {
    fn default() -> Color {
        Color {
            r: 128,
            g: 128,
            b: 128,
        }
    }
}

impl Add for Color {
    type Output = Color;

    fn add(self, rhs: Color) -> Color {
        Color {
            r: self.r.saturating_add(rhs.r),
            g: self.g.saturating_add(rhs.g),
            b: self.b.saturating_add(rhs.b),
        }
    }
}

impl Mul<f64> for Color {
    type Output = Color;

    fn mul(self, rhs: f64) -> Color {
        Color {
            r: (self.r as f64 * rhs).min(255.0) as u8,
            g: (self.g as f64 * rhs).min(255.0) as u8,
            b: (self.b as f64 * rhs).min(255.0) as u8,
        }
    }
}

impl Light {
    pub fn new() -> Light {
        Light {
            color: Color::default(),
            position: Point {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                normal: None,
            },
        }
    }

    pub fn phong(&self, point: &Point, normal: Vector, view: Vector, properties: &Properties) -> Color {
        let light_dir = Vector::from(&self.position) - Vector::from(point);
        let light_dir = light_dir.normalize();
        let ambient = self.color * properties.ambient;
        let diffuse = self.color * light_dir.dot(normal).max(0.0) * properties.diffuse;
        let reflect = (light_dir - normal * light_dir.dot(normal) * 2.).normalize();
        let specular = self.color * reflect.dot(view).max(0.0).powf(properties.shininess) * properties.specular;
        ambient + diffuse + specular
    }
}