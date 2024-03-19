extern crate bmp;

use bmp::{Image, Pixel};

use crate::arguments::AntiAliasing::SSAA;
use crate::arguments::RenderArgs;
use crate::camera::Picture;
use crate::light::SColor;

mod geometric;
mod light;
mod object;
mod camera;
mod world;
mod shapes;
mod obj_parser;
mod BVH;
mod worlds;
mod arguments;


pub fn save(picture: &Picture<SColor>, name: &str) {
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
    img.save(name).unwrap();
}


fn main() {
    // left hind system!!!
    let mut world = worlds::chess_set();
    let render_args = RenderArgs {
        height: 600,
        width: 800,
        threads: 2,
        // anti_aliasing: SSAA,
        ..RenderArgs::default()
    };
    world.render(render_args);
    save(&world.camera.picture, "output.bmp");
}
