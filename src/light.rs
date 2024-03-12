use std::ops::{Add, AddAssign, Div, Sub};
use std::ops::Mul;

use crate::geometric::{Point, Vector};
use crate::object::Properties;

#[derive(Clone, Copy)]
pub struct SColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

#[derive(Clone, Copy)]
pub struct FColor {
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

#[derive(Clone)]
pub struct Light {
    pub intensity: f64,
    pub color: SColor,
    pub position: Point,
}

pub struct VolumeLight {
    pub radiance: FColor,
    pub color: SColor,
    pub position: Point,
}

impl SColor {
    pub fn new(r: u8, g: u8, b: u8) -> SColor {
        SColor {
            r,
            g,
            b,
        }
    }
    
    pub fn black() -> SColor {
        SColor {
            r: 0,
            g: 0,
            b: 0,
        }
    }
}

impl Default for SColor {
    fn default() -> SColor {
        SColor {
            r: 0,
            g: 0,
            b: 0,
        }
    }
}

impl Add for SColor {
    type Output = SColor;

    fn add(self, rhs: SColor) -> SColor {
        SColor {
            r: self.r.saturating_add(rhs.r),
            g: self.g.saturating_add(rhs.g),
            b: self.b.saturating_add(rhs.b),
        }
    }
}

impl AddAssign for SColor {
    fn add_assign(&mut self, rhs: SColor) {
        *self = SColor {
            r: self.r.saturating_add(rhs.r),
            g: self.g.saturating_add(rhs.g),
            b: self.b.saturating_add(rhs.b),
        };
    }
}

impl Mul<f64> for SColor {
    type Output = SColor;

    fn mul(self, rhs: f64) -> SColor {
        SColor {
            r: (self.r as f64 * rhs).min(255.0) as u8,
            g: (self.g as f64 * rhs).min(255.0) as u8,
            b: (self.b as f64 * rhs).min(255.0) as u8,
        }
    }
}

impl Default for FColor {
    fn default() -> FColor {
        FColor {
            r: 0.0,
            g: 0.0,
            b: 0.0,
        }
    }
}

impl From<SColor> for Ratio {
    fn from(color: SColor) -> Ratio {
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

impl From<SColor> for FColor {
    fn from(color: SColor) -> FColor {
        FColor {
            r: color.r as f64 / 255.0,
            g: color.g as f64 / 255.0,
            b: color.b as f64 / 255.0,
        }
    }
}

impl From<FColor> for SColor {
    fn from(color: FColor) -> SColor {
        SColor {
            r: (color.r * 255.0).min(255.0) as u8,
            g: (color.g * 255.0).min(255.0) as u8,
            b: (color.b * 255.0).min(255.0) as u8,
        }
    }
}

impl Add<FColor> for FColor {
    type Output = FColor;

    fn add(self, rhs: FColor) -> FColor {
        FColor {
            r: self.r + rhs.r,
            g: self.g + rhs.g,
            b: self.b + rhs.b,
        }
    }
}

impl AddAssign<FColor> for FColor {
    fn add_assign(&mut self, rhs: FColor) {
        self.r += rhs.r;
        self.g += rhs.g;
        self.b += rhs.b;
    }
}

impl Mul<f64> for FColor {
    type Output = FColor;

    fn mul(self, rhs: f64) -> FColor {
        FColor {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        }
    }
}

impl Mul<Ratio> for FColor {
    type Output = FColor;

    fn mul(self, rhs: Ratio) -> FColor {
        FColor {
            r: self.r * rhs.r,
            g: self.g * rhs.g,
            b: self.b * rhs.b,
        }
    }
}

impl Div<f64> for FColor {
    type Output = FColor;

    fn div(self, rhs: f64) -> FColor {
        FColor {
            r: self.r / rhs,
            g: self.g / rhs,
            b: self.b / rhs,
        }
    }
}

impl FColor {
    pub fn zero() -> FColor {
        FColor {
            r: 0.0,
            g: 0.0,
            b: 0.0,
        }
    }
}


impl VolumeLight {
    pub fn new() -> VolumeLight {
        VolumeLight {
            radiance: FColor::default(),
            color: SColor::default(),
            position: Point::new(0.0, 0.0, 0.0),
        }
    }
}

impl Light {
    pub fn new() -> Light {
        Light {
            intensity: 10000.0,
            color: SColor::default(),
            position: Point::new(0.0, 0.0, 0.0),
        }
    }

    pub fn phong(&self, point: &Point, normal: Vector, view: Vector, properties: &Properties, occlusion: f64) -> SColor {
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