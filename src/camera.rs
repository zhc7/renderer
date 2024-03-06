use std::ops::{Index, IndexMut};

use crate::geometric::{Ray, Point, Vector};
use crate::light::Color;

pub struct Camera {
    pub picture: Picture,
    pub position: Point,
    pub direction: Vector,
    pub up: Vector,
    pub right: Vector,
}

impl Camera {
    pub fn point_at(&mut self, target: &Point) {
        let direction = Vector {
            x: target.x() - self.position.x(),
            y: target.y() - self.position.y(),
            z: target.z() - self.position.z(),
        };
        self.direction = direction.normalize();
        self.right = self.up.cross(self.direction).normalize();
        self.up = self.direction.cross(self.right).normalize();
    }
}

impl Camera {
    pub fn get_ray(&self, p0: u32, p1: u32) -> Ray {
        let x = (p0 as f64 + 0.5) / self.picture.width as f64;
        let y = (p1 as f64 + 0.5) / self.picture.height as f64;
        let direction = self.direction + self.right * (x - 0.5) + self.up * (0.5 - y);
        Ray {
            start: self.position.clone(),
            direction: direction.normalize(),
        }
    }
}

pub struct Picture {
    pub width: u32,
    pub height: u32,
    pub pixels: Vec<Color>,
}

impl Default for Picture {
    fn default() -> Picture {
        Picture {
            width: 800,
            height: 600,
            pixels: vec![Color::default(); 800 * 600],
        }
    }
}

impl Index<(usize, usize)> for Picture {
    type Output = Color;

    fn index(&self, index: (usize, usize)) -> &Color {
        &self.pixels[index.1 * self.width as usize + index.0]
    }
}

impl IndexMut<(usize, usize)> for Picture {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Color {
        &mut self.pixels[index.1 * self.width as usize + index.0]
    }
}

impl Camera {
    pub fn new() -> Camera {
        Camera {
            picture: Picture::default(),
            position: Point::new(0.0, 0.0, 0.0),
            direction: Vector { x: 0.0, y: 0.0, z: 1.0 },
            up: Vector { x: 0.0, y: 1.0, z: 0.0 },
            right: Vector { x: 1.0, y: 0.0, z: 0.0 },
        }
    }
}