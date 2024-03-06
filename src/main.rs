extern crate bmp;
use bmp::{Image, Pixel};
use crate::camera::Picture;
use crate::world::World;


mod geometric;
mod light;
mod object;
mod camera;
mod world;
mod shapes;


fn save(picture: &Picture) {
    let mut img = Image::new(picture.width, picture.height);
    for y in 0..picture.height {
        for x in 0..picture.width {
            let color = picture[(x as usize, y as usize)];
            img.set_pixel(x, y, Pixel {
                r: color.r,
                g: color.g,
                b: color.b,
            });
        }
    }
    img.save("output.bmp").unwrap();
}


fn main() {
    let mut world = World::new();
    world.render();
    save(&world.camera.picture);
}
