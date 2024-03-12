use std::ops::{Add, AddAssign, Div, Sub};
use std::ops::Mul;

use crate::geometric::{Point, Vector};
use crate::object::Properties;

#[derive(Clone, Copy)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

#[derive(Clone, Copy)]
pub struct Radiance {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

#[derive(Clone, Copy)]
pub struct Ratio {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

pub struct Light {
    pub radiance: Radiance,
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
            r: 0,
            g: 0,
            b: 0,
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

impl AddAssign for Color {
    fn add_assign(&mut self, rhs: Color) {
        *self = Color {
            r: self.r.saturating_add(rhs.r),
            g: self.g.saturating_add(rhs.g),
            b: self.b.saturating_add(rhs.b),
        };
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

impl Default for Radiance {
    fn default() -> Radiance {
        Radiance {
            r: 0.0,
            g: 0.0,
            b: 0.0,
        }
    }
}

impl From<Color> for Ratio {
    fn from(color: Color) -> Ratio {
        Ratio {
            r: color.r as f64 / 255.0,
            g: color.g as f64 / 255.0,
            b: color.b as f64 / 255.0,
        }
    }
}

impl Add<f64> for Ratio {
    type Output = Ratio;

    fn add(self, rhs: f64) -> Ratio {
        Ratio {
            r: self.r + rhs,
            g: self.g + rhs,
            b: self.b + rhs,
        }
    }
}

impl Mul<f64> for Ratio {
    type Output = Ratio;

    fn mul(self, rhs: f64) -> Ratio {
        Ratio {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        }
    }
}

impl Div<f64> for Ratio {
    type Output = Ratio;

    fn div(self, rhs: f64) -> Ratio {
        Ratio {
            r: self.r / rhs,
            g: self.g / rhs,
            b: self.b / rhs,
        }
    }
}

impl Add<Radiance> for Radiance {
    type Output = Radiance;

    fn add(self, rhs: Radiance) -> Radiance {
        Radiance {
            r: self.r + rhs.r,
            g: self.g + rhs.g,
            b: self.b + rhs.b,
        }
    }
}

impl AddAssign<Radiance> for Radiance {
    fn add_assign(&mut self, rhs: Radiance) {
        self.r += rhs.r;
        self.g += rhs.g;
        self.b += rhs.b;
    }
}

impl Mul<f64> for Radiance {
    type Output = Radiance;

    fn mul(self, rhs: f64) -> Radiance {
        Radiance {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        }
    }
}

impl Mul<Ratio> for Radiance {
    type Output = Radiance;

    fn mul(self, rhs: Ratio) -> Radiance {
        Radiance {
            r: self.r * rhs.r,
            g: self.g * rhs.g,
            b: self.b * rhs.b,
        }
    }
}

impl Div<f64> for Radiance {
    type Output = Radiance;

    fn div(self, rhs: f64) -> Radiance {
        Radiance {
            r: self.r / rhs,
            g: self.g / rhs,
            b: self.b / rhs,
        }
    }
}

impl Light {
    pub fn new() -> Light {
        Light {
            radiance: Radiance::default(),
            color: Color::default(),
            position: Point::new(0.0, 0.0, 0.0),
        }
    }

    pub fn phong(&self, point: &Point, normal: Vector, view: Vector, properties: &Properties, occlusion: f64) -> Color {
        assert!(1.0 - view.magnitude() < 1e-6);
        let light_dir = Vector::from(&self.position) - Vector::from(point);
        let light_dir = light_dir.normalize();
        let ambient = self.color * properties.ambient;
        if occlusion == 0.0 || light_dir.dot(normal) < 0.0 {
            return ambient;
        }
        let diffuse = self.color * light_dir.dot(normal).max(0.0) * properties.diffuse;
        let reflect = (light_dir - normal * light_dir.dot(normal) * 2.).normalize();
        let specular = self.color * reflect.dot(view).max(0.0).powf(properties.shininess) * properties.specular;
        ambient + (diffuse + specular) * occlusion
    }
}