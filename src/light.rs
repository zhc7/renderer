use std::ops::Add;
use std::ops::Mul;

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

impl Color {
    pub fn new(r: u8, g: u8, b: u8) -> Color {
        Color {
            r,
            g,
            b,
        }
    }
    
    pub fn black() -> Color {
        Color {
            r: 0,
            g: 0,
            b: 0,
        }
    }
}

impl Default for Color {
    fn default() -> Color {
        Color {
            r: 64,
            g: 64,
            b: 64,
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
            position: Point::new(0.0, 0.0, 0.0),
        }
    }

    pub fn phong(&self, point: &Point, normal: Vector, view: Vector, properties: &Properties, shadowed: bool) -> Color {
        let light_dir = Vector::from(&self.position) - Vector::from(point);
        let light_dir = light_dir.normalize();
        let ambient = self.color * properties.ambient;
        if shadowed {
            return ambient;
        }
        let diffuse = self.color * light_dir.dot(normal).max(0.0) * properties.diffuse;
        let reflect = (light_dir - normal * light_dir.dot(normal) * 2.).normalize();
        let specular = self.color * reflect.dot(view).max(0.0).powf(properties.shininess) * properties.specular;
        ambient + diffuse + specular
    }
}