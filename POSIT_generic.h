#define CV_NO_ERR 0
#define CV_NULLPTR_ERR 1
#define CV_BADSIZE_ERR 2
#define CV_OUTOFMEM_ERR 3
#define CV_BADFACTOR_ERR 4
#define CV_BADFLAG_ERR 5

#define CONVENTIONAL_LENS 6
#define FISH_EYE_LENS_STEREOGRAPHIC 7
#define FISH_EYE_LENS_EQUIDISTANCE 8
#define FISH_EYE_LENS_EQUISOLID 9
#define FISH_EYE_LENS_ORTHOGONAL 10

#include "opencv2/calib3d/calib3d_c.h"

	/* POSIT structure */
	struct POSITObject
	{
		int N;
		int ite_num;
		float* inv_matr;
		float* obj_vecs;
		float* img_vecs;
	};

    void PseudoInverse3D(float *a, float *b, int n, int method);

	int  icvCreatePOSITObject(CvPoint3D32f *points, int numPoints, POSITObject **ppObject);

	int  icvPOSIT(POSITObject *pObject, CvPoint2D32f *imagePoints, float focalLength, CvTermCriteria criteria, float xCenter, float yCenter, int lensType, float* rotation, float* translation);

	int  icvReleasePOSITObject(POSITObject ** ppObject);

	POSITObject *CreatePOSITObject(CvPoint3D32f * points, int numPoints);

	void cvPOSIT(POSITObject * pObject, CvPoint2D32f * imagePoints, double focalLength, CvTermCriteria criteria, float xCenter, float yCenter, int lensType, float* rotation, float* translation);

	void cvReleasePOSITObject(POSITObject ** ppObject);
	
	float calculate_g(float theta, int type, float g_theta);