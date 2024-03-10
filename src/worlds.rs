use crate::geometric::Point;
use crate::{light, shapes};
use crate::light::Color;
use crate::obj_parser::load;
use crate::object::Properties;
use crate::world::World;

pub fn chess_set() -> World {
    let mut world = World::new();
    // let mut object = shapes::cube();
    let mut object2 = load(r"C:\Users\bjhan\3D Objects\3D Builder\Chess set.obj".as_ref());
    // let mut object = shapes::sphere(100, 30.0);
    // let mut object2 = shapes::sphere(100, 30.0);
    let properties = Properties {
        color: Color::default(),
        ambient: 0.15,
        diffuse: 0.2,
        specular: 0.9,
        shininess: 10.0,
        refractive_index: 1.5,
        transparent: 0.8,
    };
    // object.properties = properties.clone();
    object2.properties = properties.clone();
    // world.add_object(object, Point::new(-10.0, -20.0, 120.0));
    world.add_object(object2, Point::new(0.0, -20.0, 120.0));
    world.add_light(light::Light {
        intensity: 100000.0,
        color: Color { r: 255, g: 255, b: 255 },
        position: Point::new(5.0, 200.0, 3.0),
    }, Point::new(0.0, 0.0, 0.0));
    world.camera.position = Point::new(-50.0, 60.0, 0.0);
    world.camera.point_at(&Point::new(0.0, 0.0, 120.0));
    world
}

pub fn balls() -> World {
    let mut world = World::new();
    let mut b1 = shapes::sphere(100, 10.0);
    let mut b2 = shapes::sphere(100, 15.0);
    let mut b3 = shapes::sphere(100, 25.0);
    let mut floor = shapes::cube(200.0, 200.0, 1.0);
    let properties = Properties {
        color: Color::default(),
        ambient: 0.15,
        diffuse: 0.0,
        specular: 0.9,
        shininess: 10.0,
        refractive_index: 1.5,
        transparent: 1.0,
    };
    
    b1.properties = properties.clone();
    b2.properties = properties.clone();
    b3.properties = properties.clone();
    floor.properties = Properties {
        color: Color { r: 30, g: 60, b: 200 },
        ambient: 0.15,
        diffuse: 0.8,
        specular: 0.2,
        shininess: 1.0,
        refractive_index: 1.5,
        transparent: 0.0,
    };
    world.add_object(b1, Point::new(-20.0, 0.0, 120.0));
    world.add_object(b2, Point::new(20.0, 0.0, 120.0));
    world.add_object(b3, Point::new(0.0, 0.0, 180.0));
    world.add_object(floor, Point::new(0.0, -30.0, 150.0));
    world.add_light(light::Light {
        intensity: 100000.0,
        color: Color { r: 255, g: 255, b: 255 },
        position: Point::new(100.0, 200.0, 0.0),
    }, Point::new(0.0, 0.0, 0.0));
    world.camera.position = Point::new(0.0, 30.0, 50.0);
    world.camera.point_at(&Point::new(0.0, 0.0, 120.0));
    world
}