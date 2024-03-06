use crate::camera::Camera;
use crate::light::Light;
use crate::object::Object;

pub struct Matrix {
    pub m: [[f64; 4]; 4],
}

pub struct World<'a> {
    pub objects: Vec<Object<'a>>,
    pub lights: Vec<Light>,
    pub camera: Camera,
}

impl<'a> World<'a> {
    fn new() -> World<'a> {
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
                let x = point.x * matrix.m[0][0] + point.y * matrix.m[0][1] + point.z * matrix.m[0][2] + matrix.m[0][3];
                let y = point.x * matrix.m[1][0] + point.y * matrix.m[1][1] + point.z * matrix.m[1][2] + matrix.m[1][3];
                let z = point.x * matrix.m[2][0] + point.y * matrix.m[2][1] + point.z * matrix.m[2][2] + matrix.m[2][3];
                point.x = x;
                point.y = y;
                point.z = z;
            }
        }
    }
}