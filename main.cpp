#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "SETTINGS.h"

float PI = 3.1415926535;

using namespace std;

float area(VEC3 v0, VEC3 v1){
  VEC3 cross = v0.cross(v1);
  return cross.norm()/2.0;
}

VEC3 truncate(const VEC4& v)
{
  return VEC3(v[0], v[1], v[2]);
}
VEC4 extend(const VEC3& v)
{
  return VEC4(v[0], v[1], v[2], 1.0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  fscanf(fp, "P6\n%d %d\n255\n", &xRes, &yRes);
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{

  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


class Ray{
	VEC3 ray_origin;
	VEC3 ray_direction; // NORMALIZED PLEASE
public:
	Ray(){
		ray_origin = {0,0,0};
		ray_direction = {0,0,0};
	}

	Ray(VEC3 origin, VEC3 direction){
		ray_origin = origin;
		ray_direction = direction.normalized();
	}

	VEC3 origin(){return ray_origin;}
	VEC3 direction(){return ray_direction;}

	// Returns the point on the ray for a given t
	VEC3 pointOnRay(float t){
		return ray_origin + t*ray_direction;
	}

};

class Camera{
	VEC3 cam_pos;
	VEC3 cam_gaze;
	VEC3 cam_up;
	double cam_near, cam_far, cam_fovy, cam_aspect;
	VEC3 basis_u, basis_v, basis_w;
	double screenTop, screenBottom, screenLeft, screenRight;
public:
	// Sets up camera for a 2D video (2D animation takes place at z = 0, top right is at (1920, 1080))
	Camera(int x_res, int y_res){
		cam_aspect = (float)y_res/x_res;
		cam_near = 1.0;
		cam_fovy = atan((y_res/2.0)/cam_near)*2;
		cam_pos = {x_res/2.0, y_res/2.0, 1};
		cam_gaze = {0,0,-1};
		cam_up = {0,1,0};
		cam_far = 100.0;

		// Establish camera basis
		basis_w = -cam_gaze/cam_gaze.norm();
		VEC3 upCrossW = cam_up.cross(basis_w);
		basis_u = upCrossW/upCrossW.norm();
		basis_v = basis_w.cross(basis_u);
		basis_v.normalize();

		// Screen boundaries
		screenTop = cam_near*tan(cam_fovy/2.0);
		screenBottom = -screenTop;
		screenRight = (screenTop/cam_aspect);
		screenLeft = -screenRight;

		screenTop += cam_pos[1];
		screenBottom += cam_pos[1];
		screenLeft += cam_pos[0];
		screenRight += cam_pos[0];
		//cam_near += cam_pos[2];
	}

	VEC3 position(){return cam_pos;}

	Camera(){

	}

	void printWindow(void){cout << "top = " << screenTop << ", bottom = " << screenBottom << ", right = " << screenRight << ", left = " << screenLeft << endl;}
	void printBasis(void){cout << "u = " << basis_u << ", v = " << basis_v << ", w = " << basis_w << endl;}


	Ray generateRay(int i, int j, int xRes, int yRes){
		float uCoeff = screenLeft + (screenRight-screenLeft)*((float)i+.5)/xRes;
		float vCoeff = screenBottom + (screenTop-screenBottom)*((float)j+.5)/yRes;
		VEC3 s = (vCoeff*basis_v)+(uCoeff*basis_u)-((cam_near-cam_pos[2])*basis_w);
		return Ray(cam_pos, s - cam_pos);
	}
};

class Transform_Data{
	VEC3 layer_anchor_point, layer_position, layer_scale, layer_rotation;
	float layer_opacity;
public:
	Transform_Data(){

	}

	Transform_Data(VEC3 anchor, VEC3 position, VEC3 scale, VEC3 rotation, float opacity){
		layer_anchor_point = anchor;
		layer_position = position;
		layer_scale = scale;
		layer_rotation = rotation;
		layer_opacity = opacity;
	}

	VEC3 anchor(){return layer_anchor_point;}
	VEC3 position(){return layer_position;}
	VEC3 scale(){return layer_scale;}
	VEC3 rotation(){return layer_rotation;}
	float opacity(){return layer_opacity;}

	void setAnchor(VEC3 anchor){layer_anchor_point = anchor;}
	void setPosition(VEC3 position){layer_position = position;}
	void setScale(VEC3 scale){layer_scale = scale;}
	void setRotation(VEC3 rotation){layer_rotation = rotation;}
	void setOpacity(float opacity){layer_opacity = opacity;}
};

class TransformKeyframe{
	float key_time;
	Transform_Data key_data;
public:
	TransformKeyframe(float time, Transform_Data data){
		key_time = time;
		key_data = data;
	}

	float getTime(){return key_time;}
	Transform_Data getData(){return key_data;}
};


class VidElement{
public:
	Transform_Data transform_data;
	vector<TransformKeyframe *> transform_keyframes;
	virtual float intersectRay(Ray r) = 0; // Returns t of intersection
	virtual VEC3 rayColor(Ray r) = 0; // Returns the color at the intersection of a given ray (ASSUMES THERE IS AN INTERSECTION!!! ALWAYS CHECK FIRST)
	virtual VEC3 getNormalAtIntersection(Ray r) = 0;
	virtual VidElement * applyKeyframes(float t) = 0;


	void setTransform(VEC3 anchor, VEC3 position, VEC3 scale, VEC3 rotation, float opacity){
		transform_data = Transform_Data(anchor, position, scale, rotation, opacity);
	}

	VidElement(){

	}
};


class Tri: public VidElement{
	VEC3 tri_color;
	vector<VEC3> tri_vertices;
public:
	VEC3 rayColor(Ray r){
		return tri_color;
	}

	void setAnchor(VEC3 anchor){
		transform_data.setAnchor(anchor);
	}

	void setColor(VEC3 newColor){
		tri_color = newColor;
	}

	vector<VEC3> vertices(){return tri_vertices;}
	VEC3 color(){return tri_color;}
	Transform_Data getTransformData(){return transform_data;}

	Tri(VEC3 color, vector<VEC3> vertices){
		tri_color = color;
		tri_vertices = vertices;
		// By default, the VidElement transform does nothing
		setTransform({0,0,0}, {0,0,0}, {1,1,1}, {0,0,0}, 1);
	}

	Tri(Tri * copy){
		tri_color = copy->color();
		tri_vertices = copy->vertices();
		transform_data = copy->transform_data;
	}

	VEC3 getNormalAtIntersection(Ray r){
		VEC3 v0 = tri_vertices[1] - tri_vertices[0];
		VEC3 v1 = tri_vertices[2] - tri_vertices[0];
		VEC3 n = v1.cross(v0).normalized();
		return n;
	}

	float intersectRay(Ray r){
		VEC3 n = getNormalAtIntersection(r);
		if(n.dot(r.direction()) == 0){
			return -1; // Parallel, no intersection
		}

		float t = (tri_vertices[0]-r.origin()).dot(n)/(r.direction().dot(n));
		if(t < 0){
			return -1; // Intersect behind
		}

		// Test if inside ok?
		VEC3 intersect = r.pointOnRay(t);
		VEC3 v0 = tri_vertices[1] - tri_vertices[0];
		VEC3 v1 = tri_vertices[2] - tri_vertices[0];
		VEC3 v2 = tri_vertices[2] - tri_vertices[1];
		VEC3 a = intersect - tri_vertices[0];
		VEC3 b = intersect - tri_vertices[1];
		VEC3 c = intersect - tri_vertices[2];
		float areaAll = area(v0, v1);
		float alpha = area(v2, b)/areaAll;
		float beta = area(a,v1)/areaAll;
		float gamma = area(v0, a)/areaAll;

		// Make sure the signs of cross products are the same... If not our triangle has negative area
		VEC3 areav = v0.cross(v1);
        VEC3 av = v2.cross(b);
        VEC3 bv = a.cross(v1);
        VEC3 gv = v0.cross(a);
        float asign = av[2]*areav[2];
        float bsign = bv[2]*areav[2];
        float gsign = gv[2]*areav[2];

		if((asign >= 0.0) && (bsign >= 0.0) && (gsign >= 0.0) && (alpha >= 0) && (alpha <= 1) && (beta >= 0) && (beta <= 1) && (gamma >= 0) && (gamma <= 1)){
			return t;
		}else{
			return -1;
		}
	}

	void rotate(VEC3 rotate){
		float rad_x = rotate[0];
		float rad_y = rotate[1];
		float rad_z = rotate[2];
		MATRIX4 RX;
		RX.setZero();
		RX(0,0) = 1;
		RX(1,1) = cos(rad_x);
		RX(1,2) = -sin(rad_x);
		RX(2,1) = sin(rad_x);
		RX(2,2) = cos(rad_x);
		RX(3,3) = 1;

		MATRIX4 RY;
		RY.setZero();
		RY(1,1) = 1;
		RY(0,0) = cos(rad_y);
		RY(2,0) = -sin(rad_y);
		RY(0,2) = sin(rad_y);
		RY(2,2) = cos(rad_y);
		RY(3,3) = 1;

		MATRIX4 RZ;
		RZ.setZero();
		RZ(2,2) = 1;
		RZ(0,0) = cos(rad_z);
		RZ(0,1) = -sin(rad_z);
		RZ(1,0) = sin(rad_z);
		RZ(1,1) = cos(rad_z);
		RZ(3,3) = 1;

		MATRIX4 T1;
		T1.setZero();
		T1(0,0) = 1;
		T1(1,1) = 1;
		T1(2,2) = 1;
		T1(3,3) = 1;
		T1(0,3) = -transform_data.anchor()[0];
		T1(1,3) = -transform_data.anchor()[1];
		T1(2,3) = -transform_data.anchor()[2];

		MATRIX4 T2;
		T2.setZero();
		T2(0,0) = 1;
		T2(1,1) = 1;
		T2(2,2) = 1;
		T2(3,3) = 1;
		T2(0,3) = transform_data.anchor()[0];
		T2(1,3) = transform_data.anchor()[1];
		T2(2,3) = transform_data.anchor()[2];

		vector<VEC3> newVertices;
		for(int i = 0; i < 3; i++){
			newVertices.push_back(truncate(T2*RX*RY*RZ*T1*extend(tri_vertices[i])));
		}
		tri_vertices = newVertices;
	}

	void translateAnchor(VEC3 newAnchor){
		transform_data.setAnchor(newAnchor);
	}

	void translate(VEC3 translate){
		float trans_x = translate[0];
		float trans_y = translate[1];
		float trans_z = translate[2];
		MATRIX4 T;
		T.setZero();
		T(0,0) = 1;
		T(1,1) = 1;
		T(2,2) = 1;
		T(3,3) = 1;
		T(0,3) = trans_x;
		T(1,3) = trans_y;
		T(2,3) = trans_z;

		vector<VEC3> newVertices;
		for(int i = 0; i < 3; i++){
			newVertices.push_back(truncate(T*extend(tri_vertices[i])));
		}
		tri_vertices = newVertices;
	}

	void scale(VEC3 scale){
		float scale_x = scale[0];
		float scale_y = scale[1];
		float scale_z = scale[2];
		MATRIX4 T1;
		T1.setZero();
		T1(0,0) = 1;
		T1(1,1) = 1;
		T1(2,2) = 1;
		T1(3,3) = 1;
		T1(0,3) = -transform_data.anchor()[0];
		T1(1,3) = -transform_data.anchor()[1];
		T1(2,3) = -transform_data.anchor()[2];

		MATRIX4 T2;
		T2.setZero();
		T2(0,0) = 1;
		T2(1,1) = 1;
		T2(2,2) = 1;
		T2(3,3) = 1;
		T2(0,3) = transform_data.anchor()[0];
		T2(1,3) = transform_data.anchor()[1];
		T2(2,3) = transform_data.anchor()[2];

		MATRIX4 SCALE;
		SCALE.setZero();
		SCALE(0,0) = scale_x;
		SCALE(1,1) = scale_y;
		SCALE(2,2) = scale_z;
		SCALE(3,3) = 1;

		vector<VEC3> newVertices;
		for(int i = 0; i < 3; i++){
			newVertices.push_back(truncate(T2*SCALE*T1*extend(tri_vertices[i])));
		}
		tri_vertices = newVertices;
	}

	Tri * translated(VEC3 translate){
		Tri * out = new Tri(this);
		out->translate(translate);
		return out;
	}

	void newKeyFrame(float time, Transform_Data data){
		transform_keyframes.push_back(new TransformKeyframe(time, data));
	}

	// Returns a new VidElement with transformations according to t
	VidElement * applyKeyframes(float t){
		// Assumes keyframes are in order
		cout << "t = " << t;
		TransformKeyframe * last_keyframe;
		TransformKeyframe * next_keyframe;
		int i = 0;
		while((i < transform_keyframes.size()) && (t > transform_keyframes[i]->getTime())){
			i++;
		}


		float percent = 1;
		if(i == 0){
			last_keyframe = transform_keyframes[i]; // Haven't reached the first keyframe yet, keep percent = 1
			next_keyframe = transform_keyframes[i];
		}else if(i == transform_keyframes.size()){
			last_keyframe = transform_keyframes[i-1]; // Already past the last keyframe, keep percent = 1
			next_keyframe = transform_keyframes[i-1];
		}else{
			last_keyframe = transform_keyframes[i-1];
			next_keyframe = transform_keyframes[i];
			// Just do linear interpolation for now
			float total_dt = next_keyframe->getTime() - last_keyframe->getTime();
			float current_dt = t - last_keyframe->getTime();
			percent = current_dt/total_dt;
		}

		cout << " last time = " << last_keyframe->getTime() << ", next time = " << next_keyframe->getTime();
		cout << "percent = " << percent << endl;

		Transform_Data last_data = last_keyframe->getData();
		Transform_Data next_data = next_keyframe->getData();

		VEC3 d_pos = next_data.position() - last_data.position();
		VEC3 d_rot = next_data.rotation() - last_data.rotation();
		VEC3 d_scale = next_data.scale() - last_data.scale();
		VEC3 d_anchor = next_data.anchor() - last_data.anchor();
		float d_opacity = next_data.opacity() - last_data.opacity();

		Tri * out = new Tri(this);
		out->translateAnchor(last_data.anchor() + percent*d_anchor);
		out->translate(last_data.position() + percent*d_pos);
		out->rotate(last_data.rotation() + percent*d_rot);
		out->scale(last_data.scale() + percent*d_scale);

		return out;
	}

};

class Scene{
	Camera scene_cam;
	vector<VidElement*> scene_elements_original;
	vector<VidElement*> scene_elements;
public:
	Scene(Camera cam, vector<VidElement*> elements){
		scene_cam = cam;
		scene_elements_original = elements;
	}

	Scene(){

	}

	VidElement * closestElement(Ray r){
		double closestDist = -1.0;
		VidElement * closestElement = NULL;

		for(VidElement * v : scene_elements){
			double dist = v->intersectRay(r);
			if(dist > 0){
				if((closestDist < 0) || (dist < closestDist)){
					closestDist = dist;
					closestElement = v;
				}
			}
		}

		return closestElement;
	}


	VEC3 rayColor(Ray r){
		VidElement * closest = closestElement(r);
		if(closest){
			return closest->rayColor(r);
		}else{
			return {0,0,0};
		}
	}

	void applyKeyframes(float t){
		scene_elements.clear();
		for(VidElement * v : scene_elements_original){
			scene_elements.push_back(v->applyKeyframes(t));
		}
	}

	float * rasterize(int xRes, int yRes, double t){
		applyKeyframes(t);
		float * pixels = new float[3*xRes*yRes];
		Ray r;
		for(int y = 0; y < yRes; y++){
			for(int x = 0; x < xRes; x++){
				int index = ((yRes-1-y)*xRes + x);
				r = scene_cam.generateRay(x, y, xRes, yRes);
				VEC3 color = rayColor(r);
				pixels[(index*3)+0] = color[0]*255;
				pixels[(index*3)+1] = color[1]*255;
				pixels[(index*3)+2] = color[2]*255;
			}
		}
		return pixels;
	}
};

class Renderer{
	float fps;
	int xRes, yRes;
	Scene scene;
	string path;
public:
	Renderer(string path, Scene scene, int xRes, int yRes, float fps){
		this->fps = fps;
		this->xRes = xRes;
		this->yRes = yRes;
		this->scene = scene;
		this->path = path;
	}

	void render(int start_frame, int end_frame){
		int digits = 0;
		int temp = end_frame;
		while(temp > 0){
			temp /= 10;
			digits++;
		}

		for(int i = start_frame; i < end_frame; i++){
			if((i%10) == 0){
				//cout << i << endl;
			}
			float * pixels = scene.rasterize(xRes, yRes, i*(1.0/fps));
			stringstream ss;
			ss << setw(digits) << setfill('0') << i << ".ppm";
			string filename = path+ss.str();
			writePPM(filename, xRes, yRes, pixels);
		}
	}

};


int main(int argc, char** argv)
{	
	int xRes = 1920;
	int yRes = 1080;

	Camera camera = Camera(xRes, yRes);
	vector<VidElement*> elements;

	Tri * tri1 = new Tri({.25,.25,1}, {{0,1080,0}, {0,0,0}, {1920,0,0}});
	tri1->setAnchor({1920/2.0, 1080/2.0, 0});
	tri1->newKeyFrame(.5, Transform_Data({0,0,0}, {1920/2.0, 1080/2.0, 0}, {1,1,1}, {0,0,0}, 1));
	tri1->newKeyFrame(1, Transform_Data({0,0,0}, {1920/2.0, 1080/2.0, 0}, {1,1,1}, {PI/4,0,0}, 1));

	elements.push_back(tri1);
	Scene scene = Scene(camera, elements);
	Renderer r = Renderer("oven/", scene, xRes, yRes, 30);
	r.render(0,35);

  	return 0;
}
