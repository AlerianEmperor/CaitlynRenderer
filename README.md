<img width="1920" height="1080" alt="11 9 fps" src="https://github.com/user-attachments/assets/4f94f6be-bc10-4895-9c62-e4bca468460b" /># Caitlyn Renderer
A Real Time Path Tracer in GLSL.

![rank3police](https://github.com/AlerianEmperor/CaitlynRenderer/assets/93391908/aa5cb991-f7e8-4580-9e90-1ef116486bd6)

![DiningRoom](https://github.com/AlerianEmperor/CaitlynRenderer/assets/93391908/6992d4d8-f3e8-49da-92aa-c0be57f0d19b)

![iron](https://github.com/AlerianEmperor/CaitlynRenderer/assets/93391908/dbf68cb0-7664-4114-be96-e84a4a20db7f)

![cwbvh_1024](https://github.com/AlerianEmperor/CaitlynRenderer/assets/93391908/f0617806-9d95-4353-9f7c-41e94c25aec8)

These images were rendered by my 10 years old laptop, if I have a new laptop or GPU, the performance will be much higher!

Final update:
convert to compute shader, but it is not so much faster than expected, wbvh is a hoax also, may review cwbvh in the future. 

This will be the last update before realsing my newest game engine.

over all Cailtyn is now 3.3 time faster than my old version.

Ray tracing is so boring if you do it for a while, game engine is so much better!

fireroom 144k triangles, 1024 x 720, 3 depth, nee, 11.9 fps, sbvh

<img width="1920" height="1080" alt="11 9 fps" src="https://github.com/user-attachments/assets/16b0930c-408d-4599-9310-85f16a87f1da" />

# Control
W, A, S, D to move, left mouse to navigate, right mouse to adjust UI.

# Features
- Real time path tracing.
- Unidirectional Path Tracing.
- SAH BVH.
- SBVH (Spatial Split BVH), which run 20% - 30% faster than SAH BVH.
- CWBVH (Compressed Wide BVH), which run 2 to 4 times faster than SBVH.
- Disney BSDF.
- Interactive control and update scene setting using GUI.

<i
