use std::hash::Hash;
use std::ops::{Add, Div, Index, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub normal: Option<Vector>,
}

#[derive(Debug, Clone, Hash)]
pub struct Segment {
    pub start: Point,
    pub end: Point,
}

pub struct Ray {
    pub start: Point,
    pub direction: Vector,
}

#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub struct Triangle<'a> {
    pub a: &'a Point,
    pub b: &'a Point,
    pub c: &'a Point,
    pub normal: Option<Vector>,
}

impl Vector {
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Vector {
        let magnitude = self.magnitude();
        Vector {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
        }
    }

    pub fn dot(&self, rhs: Vector) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn cross(&self, rhs: Vector) -> Vector {
        Vector {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }
}

impl PartialEq<Self> for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

impl Eq for Vector {}


impl Hash for Vector {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.to_bits().hash(state);
        self.y.to_bits().hash(state);
        self.z.to_bits().hash(state);
    }
}

impl From<&Point> for Vector {
    fn from(point: &Point) -> Vector {
        Vector {
            x: point.x,
            y: point.y,
            z: point.z,
        }
    }
}

impl Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Vector {
        Vector {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Vector {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Vector {
        Vector {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, rhs: f64) -> Vector {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for Vector {
    type Output = Vector;

    fn div(self, rhs: f64) -> Vector {
        Vector {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Point {
            x,
            y,
            z,
            normal: None,
        }
    }
}

impl PartialEq<Self> for Point {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

impl Eq for Point {}

impl Hash for Point {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.to_bits().hash(state);
        self.y.to_bits().hash(state);
        self.z.to_bits().hash(state);
    }
}

impl From<Vector> for Point {
    fn from(vector: Vector) -> Point {
        Point {
            x: vector.x,
            y: vector.y,
            z: vector.z,
            normal: None,
        }
    }
}

impl<'a> Triangle<'a> {
    pub fn new(a: &'a Point, b: &'a Point, c: &'a Point) -> Triangle<'a> {
        Triangle {
            a,
            b,
            c,
            normal: None,
        }
    }

    pub fn calc_norm(&mut self) {
        let ab = Vector::from(self.b) - Vector::from(self.a);
        let ac = Vector::from(self.c) - Vector::from(self.a);
        let normal = ab.cross(ac).normalize();
        self.normal = Some(normal);
    }

    pub fn size(&self) -> f64 {
        let ab = Vector::from(self.b) - Vector::from(self.a);
        let ac = Vector::from(self.c) - Vector::from(self.a);
        let cross = ab.cross(ac);
        cross.magnitude() / 2.0
    }

    pub fn intersect(&self, ray: &Ray) -> Option<Point> {
        let ab = Vector::from(self.b) - Vector::from(self.a);
        let ac = Vector::from(self.c) - Vector::from(self.a);
        let normal = ab.cross(ac).normalize();
        let d = -normal.dot(Vector::from(self.a));
        let t = -(normal.dot(Vector::from(&ray.start)) + d) / normal.dot(ray.direction);
        if t < 0.0 {
            return None;
        }
        let point = Vector::from(&ray.start) + ray.direction * t;
        let ap = point - Vector::from(self.a);
        let dot00 = ac.dot(ac);
        let dot01 = ac.dot(ab);
        let dot02 = ac.dot(ap);
        let dot11 = ab.dot(ab);
        let dot12 = ab.dot(ap);
        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
            Some(Point::from(point))
        } else {
            None
        }
    }
}

impl<'a> Index<usize> for Triangle<'a> {
    type Output = Point;

    fn index(&self, index: usize) -> &Point {
        match index {
            0 => &self.a,
            1 => &self.b,
            2 => &self.c,
            _ => panic!("Index out of bounds"),
        }
    }
}

