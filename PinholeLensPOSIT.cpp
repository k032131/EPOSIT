#include "POSIT_pinhole.h"
#include "stdio.h"
#include <iostream>

using namespace std;
using namespace cv;

int main()
{
	//POSIT算法不是以图像左上角为圆点的，而是以u0,v0为圆点，而且x指的column，y指的Row
	//pinhole lens
	float x_Point0 = 876.3492 - 807.41661;
	float y_Point0 = 671.5425 - 602.4737;
	float x_Point1 = 873.7555 - 807.41661;
	float y_Point1 = 673.3157 - 602.4737;
	float x_Point2 = 905.1681 - 807.41661;
	float y_Point2 = 685.1693 - 602.4737;
	float x_Point3 = 865.9635 - 807.41661;
	float y_Point3 = 704.5122 - 602.4737;

	
	double focal_Length = 6 / 5.3 * 1200;//7.18 * 1600;////
	CvPOSITObject* positObject;
	std::vector<CvPoint3D32f> modelPoints;
	modelPoints.push_back(cvPoint3D32f(0.0f, 0.0f, 0.0f)); //The first must be (0,0,0)
	modelPoints.push_back(cvPoint3D32f(0.0f, 0.0f, -10.0f));
	modelPoints.push_back(cvPoint3D32f(50.0f, 0.0f, -10.0f));
	modelPoints.push_back(cvPoint3D32f(0.0f, 50.0f, -10.0f));
	

	positObject = cvCreatePOSITObject(&modelPoints[0], (int)modelPoints.size());

	std::vector<CvPoint2D32f> srcImagePoints;
	srcImagePoints.push_back(cvPoint2D32f(x_Point0, y_Point0));
	srcImagePoints.push_back(cvPoint2D32f(x_Point1, y_Point1));
	srcImagePoints.push_back(cvPoint2D32f(x_Point2, y_Point2));
	srcImagePoints.push_back(cvPoint2D32f(x_Point3, y_Point3));
	

	float* rotation_matrix = new float[9];
	float* translation_vector = new float[3];
	CvTermCriteria criteria = cvTermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 50000, 1.0e-4f);
	cvPOSIT(positObject, &srcImagePoints[0], focal_Length, criteria, rotation_matrix, translation_vector);
	std::cout << "Rotation Matrix:\n" << rotation_matrix[0] << "," << rotation_matrix[1] << "," << rotation_matrix[2] << "\n" << rotation_matrix[3] \
		<< "," << rotation_matrix[4] << "," << rotation_matrix[5] << "\n" << rotation_matrix[6] << "," << rotation_matrix[7] << "," << rotation_matrix[8] << "\n";
	std::cout << "Rotate Angle:" << atan2(rotation_matrix[7], rotation_matrix[8]) * 180 / 3.14 << "," << atan2(-rotation_matrix[6], sqrt(rotation_matrix[7] * rotation_matrix[7] + rotation_matrix[8] * rotation_matrix[8])) * 180 / 3.14 << "," << atan2(rotation_matrix[3], rotation_matrix[0]) * 180 / 3.14 << "\n";
	std::cout << "Translation Matrix:" << translation_vector[0] << "," << translation_vector[1] << "," << translation_vector[2] << "\n";

	system("pause");
	return CV_NO_ERR;
}

