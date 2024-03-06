use std::ops::{Index, IndexMut};
use crate::geometric::{Vector, Ray, Point};
use crate::light::{Color, Light};
use crate::object::Object;

pub struct Camera {
    position: Point,
    direction: Vector,
}

struct Picture {
    width: u32,
    height: u32,
    pixels: Vec<Color>,
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
            position: Point::new(0.0, 0.0, 0.0),
            direction: Vector { x: 0.0, y: 0.0, z: 1.0 },
        }
    }
}