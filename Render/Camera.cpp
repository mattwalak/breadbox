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

		if((i == 960) && (j == 540)){
			cout << "origin = " << cam_pos << "s = " << s << endl;
		}
		return Ray(cam_pos, s - cam_pos);
	}
};