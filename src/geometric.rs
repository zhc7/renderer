use std::cell::RefCell;
use std::hash::Hash;
use std::ops::{Add, Div, Index, Mul, Neg, Sub};
use std::rc::Rc;

use crate::world::Matrix;

#[derive(Debug, Clone, Copy)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector {
            x,
            y,
            z,
        }
    }

    pub fn zero() -> Vector {
        Vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
struct _Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub normal: Option<Vector>,
}

#[derive(Debug)]
pub struct Point {
    inner : Rc<RefCell<_Point>>,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Point {
            inner: Rc::new(RefCell::new(_Point {
                x,
                y,
                z,
                normal: None,
            })),
        }
    }
    
    pub fn origin() -> Point {
        Point::new(0., 0., 0.)
    }

    pub fn x(&self) -> f64 {
        self.inner.borrow().x
    }

    pub fn y(&self) -> f64 {
        self.inner.borrow().y
    }

    pub fn z(&self) -> f64 {
        self.inner.borrow().z
    }

    pub fn normal(&self) -> Option<Vector> {
        self.inner.borrow().normal
    }

    pub fn set_normal(&self, normal: Vector) {
        self.inner.borrow_mut().normal = Some(normal);
    }

    pub fn transform(&self, matrix: &Matrix) {
        self.inner.borrow_mut().transform(matrix);
    }

    pub fn distance(&self, other: &Point) -> f64 {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        let dz = self.z() - other.z();
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl Clone for Point {
    fn clone(&self) -> Point {
        Point {
            inner: Rc::clone(&self.inner),
        }
    }
}

impl Hash for Point {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.inner.borrow().hash(state)
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.inner.borrow().eq(&*other.inner.borrow())
    }
}

impl Eq for Point {}

#[derive(Debug, Clone, Hash)]
pub struct Segment {
    pub start: Point,
    pub end: Point,
}

pub struct Ray {
    pub start: Point,
    pub direction: Vector,
}

#[derive(Debug, Clone, PartialEq, Hash)]
pub struct Triangle {
    pub a: Point,
    pub b: Point,
    pub c: Point,
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

    pub fn transform(&mut self, matrix: &Matrix) {
        let x = self.x * matrix.m[0][0] + self.y * matrix.m[1][0] + self.z * matrix.m[2][0];
        let y = self.x * matrix.m[0][1] + self.y * matrix.m[1][1] + self.z * matrix.m[2][1];
        let z = self.x * matrix.m[0][2] + self.y * matrix.m[1][2] + self.z * matrix.m[2][2];
        self.x = x;
        self.y = y;
        self.z = z;
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

impl From<&_Point> for Vector {
    fn from(point: &_Point) -> Vector {
        Vector {
            x: point.x,
            y: point.y,
            z: point.z,
        }
    }
}

impl From<&Point> for Vector {
    fn from(point: &Point) -> Vector {
        Vector {
            x: point.x(),
            y: point.y(),
            z: point.z(),
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

impl _Point {
    pub fn new(x: f64, y: f64, z: f64) -> _Point {
        _Point {
            x,
            y,
            z,
            normal: None,
        }
    }

    pub fn transform(&mut self, matrix: &Matrix) {
        let x = self.x * matrix.m[0][0] + self.y * matrix.m[0][1] + self.z * matrix.m[0][2] + matrix.m[0][3];
        let y = self.x * matrix.m[1][0] + self.y * matrix.m[1][1] + self.z * matrix.m[1][2] + matrix.m[1][3];
        let z = self.x * matrix.m[2][0] + self.y * matrix.m[2][1] + self.z * matrix.m[2][2] + matrix.m[2][3];
        self.x = x;
        self.y = y;
        self.z = z;
    }

    pub fn distance(&self, other: &_Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl PartialEq<Self> for _Point {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

impl Eq for _Point {}

impl Hash for _Point {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.to_bits().hash(state);
        self.y.to_bits().hash(state);
        self.z.to_bits().hash(state);
    }
}

impl From<Vector> for _Point {
    fn from(vector: Vector) -> _Point {
        _Point {
            x: vector.x,
            y: vector.y,
            z: vector.z,
            normal: None,
        }
    }
}

impl From<Vector> for Point {
    fn from(vector: Vector) -> Point {
        Point::new(vector.x, vector.y, vector.z)
    }
}

impl Triangle {
    pub fn new(a: &Point, b: &Point, c: &Point) -> Triangle {
        Triangle {
            a: a.clone(),
            b: b.clone(),
            c: c.clone(),
            normal: None,
        }
    }

    pub fn calc_norm(&mut self) {
        let ab = Vector::from(&self.b) - Vector::from(&self.a);
        let ac = Vector::from(&self.c) - Vector::from(&self.a);
        let normal = ab.cross(ac).normalize();
        self.normal = Some(normal);
    }

    pub fn size(&self) -> f64 {
        let ab = Vector::from(&self.b) - Vector::from(&self.a);
        let ac = Vector::from(&self.c) - Vector::from(&self.a);
        let cross = ab.cross(ac);
        cross.magnitude() / 2.0
    }

    pub fn intersect(&self, ray: &Ray) -> Option<(f64, Point)> {
        let ab = Vector::from(&self.b) - Vector::from(&self.a);
        let ac = Vector::from(&self.c) - Vector::from(&self.a);
        let normal = self.normal.unwrap();
        let d = -normal.dot(Vector::from(&self.a));
        let t = -(normal.dot(Vector::from(&ray.start)) + d) / normal.dot(ray.direction);
        if t < 0.0 {
            return None;
        }
        let point = Vector::from(&ray.start) + ray.direction * t;
        let ap = point - Vector::from(&self.a);
        let dot00 = ac.dot(ac);
        let dot01 = ac.dot(ab);
        let dot02 = ac.dot(ap);
        let dot11 = ab.dot(ab);
        let dot12 = ab.dot(ap);
        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
            Some((t, Point::from(point)))
        } else {
            None
        }
    }
}

impl Index<usize> for Triangle {
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

