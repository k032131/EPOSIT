#include "POSIT_generic.h"

/*F///////////////////////////////////////////////////////////////////////////////////////
//    Name:       icvPseudoInverse3D
//    Purpose:    Pseudoinverse N x 3 matrix     N >= 3
//    Context:
//    Parameters:
//                a - input matrix
//                b - pseudoinversed a
//                n - number of rows in a
//                method - if 0, then b = inv(transpose(a)*a) * transpose(a)
//                         if 1, then SVD used.
//    Returns:
//    Notes:      Both matrix are stored by n-dimensional vectors.
//                Now only method == 0 supported.
//F*/
void PseudoInverse3D(float *a, float *b, int n, int method)
{
	if (method == 0)
	{
		float ata00 = 0;
		float ata11 = 0;
		float ata22 = 0;
		float ata01 = 0;
		float ata02 = 0;
		float ata12 = 0;

		int k;
		/* compute matrix ata = transpose(a) * a  */
		for (k = 0; k < n; k++)
		{
			float a0 = a[k];
			float a1 = a[n + k];
			float a2 = a[2 * n + k];

			ata00 += a0 * a0;
			ata11 += a1 * a1;
			ata22 += a2 * a2;

			ata01 += a0 * a1;
			ata02 += a0 * a2;
			ata12 += a1 * a2;
		}
		/* inverse matrix ata */
		{
			float p00 = ata11 * ata22 - ata12 * ata12;
			float p01 = -(ata01 * ata22 - ata12 * ata02);
			float p02 = ata12 * ata01 - ata11 * ata02;

			float p11 = ata00 * ata22 - ata02 * ata02;
			float p12 = -(ata00 * ata12 - ata01 * ata02);
			float p22 = ata00 * ata11 - ata01 * ata01;

			float det = 0;
			det += ata00 * p00;
			det += ata01 * p01;
			det += ata02 * p02;

			const float inv_det = 1 / det;

			/* compute resultant matrix */
			for (k = 0; k < n; k++)
			{
				float a0 = a[k];
				float a1 = a[n + k];
				float a2 = a[2 * n + k];

				b[k] = (p00 * a0 + p01 * a1 + p02 * a2) * inv_det;
				b[n + k] = (p01 * a0 + p11 * a1 + p12 * a2) * inv_det;
				b[2 * n + k] = (p02 * a0 + p12 * a1 + p22 * a2) * inv_det;
			}
		}
	}

	/*if ( method == 1 )
	{
	}
	*/

	return;
}

float calculate_g(float theta, int type, float g_theta)
{
	switch (type)
	{
	case 6:
	{
		g_theta = cos(theta);
		return g_theta;
		break;
	}
	case 7:
	{
		g_theta = cos(theta / 2)*cos(theta / 2);
		return g_theta;
		break;
	}

	case 8:
	{
		g_theta = sin(theta) / theta;
		return g_theta;
		break;
	}


	case 9:
	{
		g_theta = cos(theta / 2);
		return g_theta;
		break;
	}

	case 10:
	{
		g_theta = 1.0;
		return g_theta;
		break;
	}
	}
}

int  icvCreatePOSITObject(CvPoint3D32f *points, int numPoints,  POSITObject **ppObject)
{
	int i;

	/* Compute size of required memory */
	/* buffer for inverse matrix = N*3*float */
	/* buffer for storing weakImagePoints = numPoints * 2 * float */
	/* buffer for storing object vectors = N*3*float */
	/* buffer for storing image vectors = N*2*float */

	int N = numPoints - 1;
	int inv_matr_size = N * 3 * sizeof(float);
	int obj_vec_size = inv_matr_size;
	int img_vec_size = N * 2 * sizeof(float);
	POSITObject *pObject;

	/* check bad arguments */
	if (points == NULL)
		return CV_NULLPTR_ERR;
	if (numPoints < 4)
		return CV_BADSIZE_ERR;
	if (ppObject == NULL)
		return CV_NULLPTR_ERR;

	/* memory allocation */
	pObject = (POSITObject *)cvAlloc(sizeof(POSITObject) +
		inv_matr_size + obj_vec_size + img_vec_size);

	if (!pObject)
		return CV_OUTOFMEM_ERR;

	/* part the memory between all structures */
	pObject->N = N;
	pObject->ite_num = 0;
	pObject->inv_matr = (float *)((char *)pObject + sizeof(POSITObject));
	pObject->obj_vecs = (float *)((char *)(pObject->inv_matr) + inv_matr_size);
	pObject->img_vecs = (float *)((char *)(pObject->obj_vecs) + obj_vec_size);

	/****************************************************************************************\
	*          Construct object vectors from object points     3*N matrix                              *
	\****************************************************************************************/
	for (i = 0; i < numPoints - 1; i++)
	{
		pObject->obj_vecs[i] = points[i + 1].x - points[0].x;
		pObject->obj_vecs[N + i] = points[i + 1].y - points[0].y;
		pObject->obj_vecs[2 * N + i] = points[i + 1].z - points[0].z;
	}
	/****************************************************************************************\
	*   Compute pseudoinverse matrix                                                         *
	\****************************************************************************************/
	PseudoInverse3D(pObject->obj_vecs, pObject->inv_matr, N, 0);

	*ppObject = pObject;
	return CV_NO_ERR;
}


int  icvPOSIT(POSITObject *pObject, CvPoint2D32f *imagePoints, float focalLength, CvTermCriteria criteria, float xCenter, float yCenter, int lensType, float* rotation, float* translation)
{
	int i, j, k;
	int count = 0, converged = 0;
	float scale = 0, inv_Z = 0;
	float diff = (float)criteria.epsilon;

	/* Check bad arguments */
	if (imagePoints == NULL)
		return CV_NULLPTR_ERR;
	if (pObject == NULL)
		return CV_NULLPTR_ERR;
	if (focalLength <= 0)
		return CV_BADFACTOR_ERR;
	if (!rotation)
		return CV_NULLPTR_ERR;
	if (!translation)
		return CV_NULLPTR_ERR;
	if ((criteria.type == 0) || (criteria.type > (CV_TERMCRIT_ITER | CV_TERMCRIT_EPS)))
		return CV_BADFLAG_ERR;
	if ((criteria.type & CV_TERMCRIT_EPS) && criteria.epsilon < 0)
		return CV_BADFACTOR_ERR;
	if ((criteria.type & CV_TERMCRIT_ITER) && criteria.max_iter <= 0)
		return CV_BADFACTOR_ERR;

	/* init variables */
	float inv_focalLength = 1 / focalLength;
	int N = pObject->N;
	float *objectVectors = pObject->obj_vecs;
	float *invMatrix = pObject->inv_matr;
	float *imgVectors = pObject->img_vecs;
	float _g_M0 = 0;
	float _g_Mi = 0;
	float g_M0 = 0;
	float g_Mi = 0;
	
	while (!converged)
	{   
		float theta_M0 = asin(sqrt((imagePoints[0].x - xCenter)*(imagePoints[0].x - xCenter) + \
			(imagePoints[0].y - yCenter)*(imagePoints[0].y - yCenter)) / focalLength);
		g_M0 = calculate_g(theta_M0, lensType, _g_M0);
		if (count == 0)
		{
			/* subtract out origin to get image vectors */
			for (i = 0; i < N; i++)
			{
				float theta_Mi = asin(sqrt((imagePoints[i+1].x - xCenter)*(imagePoints[i+1].x - xCenter) + \
					(imagePoints[i + 1].y - yCenter)*(imagePoints[i + 1].y - yCenter)) / focalLength);
				g_Mi = calculate_g(theta_Mi, lensType, _g_Mi);
				imgVectors[i] = imagePoints[i + 1].x / cos(theta_Mi)* g_Mi - imagePoints[0].x / cos(theta_M0) * g_M0;//xi-x0
				imgVectors[N + i] = imagePoints[i + 1].y / cos(theta_Mi) * g_Mi - imagePoints[0].y / cos(theta_M0) * g_M0;//yi-y0
			}
		}
		else
		{
			diff = 0;//即|ibuxinuo(n)-ibuxinuo(n-1)|
			/* Compute new SOP (scaled orthograthic projection) image from pose */
			for (i = 0; i < N; i++)
			{
				float theta_Mi = asin(sqrt((imagePoints[i + 1].x - xCenter)*(imagePoints[i + 1].x - xCenter) + \
					(imagePoints[i + 1].y - yCenter)*(imagePoints[i + 1].y - yCenter)) / focalLength);
				g_Mi = calculate_g(theta_Mi, lensType, _g_Mi);
				/* objectVector * k */ //即M0Mi.k
				float old;
				//计算M0Mi.k，其中objectVectors[i]/[N+i]/[2*N+i]分别代表世界坐标中的向量x,y,z坐标，rotation[6]/[7]/[8]代表坐标轴k
				float tmp = objectVectors[i] * rotation[6] /*[2][0]*/ +
					objectVectors[N + i] * rotation[7]     /*[2][1]*/ +
					objectVectors[2 * N + i] * rotation[8] /*[2][2]*/;

				tmp *= inv_Z;//计算ibuxinuo(i)  = 1/Z0 * M0Mi.k
				tmp += 1;//1+ibuxinuo(i)

				old = imgVectors[i];//上一次的x'
				imgVectors[i] = (imagePoints[i + 1].x - xCenter) * tmp / cos(theta_Mi) * g_Mi - (imagePoints[0].x - xCenter) / cos(theta_M0) * g_M0;//x'

				diff = MAX(diff, (float)fabs(imgVectors[i] - old));

				old = imgVectors[N + i];//上一次的y'
				imgVectors[N + i] = (imagePoints[i + 1].y - yCenter) * tmp / cos(theta_Mi) * g_Mi - (imagePoints[0].y - yCenter) / cos(theta_M0) * g_M0;//y'

				diff = MAX(diff, (float)fabs(imgVectors[N + i] - old));
			}
		}

		/* calculate I and J vectors */
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 3; j++)
			{
				rotation[3 * i + j] /*[i][j]*/ = 0;
				for (k = 0; k < N; k++)
				{
					rotation[3 * i + j] /*[i][j]*/ += invMatrix[j * N + k] * imgVectors[i * N + k];
				}
			}
		}

		float inorm =//I
			rotation[0] /*[0][0]*/ * rotation[0] /*[0][0]*/ +
			rotation[1] /*[0][1]*/ * rotation[1] /*[0][1]*/ +
			rotation[2] /*[0][2]*/ * rotation[2] /*[0][2]*/;

		float jnorm =//J
			rotation[3] /*[1][0]*/ * rotation[3] /*[1][0]*/ +
			rotation[4] /*[1][1]*/ * rotation[4] /*[1][1]*/ +
			rotation[5] /*[1][2]*/ * rotation[5] /*[1][2]*/;

		const float invInorm = cvInvSqrt(inorm) / cos(theta_M0) * g_M0;//平方根的倒数，即1/s1
		const float invJnorm = cvInvSqrt(jnorm) / cos(theta_M0) * g_M0;//1/s2

		inorm *= invInorm;//i
		jnorm *= invJnorm;//j

		rotation[0] /*[0][0]*/ *= invInorm;
		rotation[1] /*[0][1]*/ *= invInorm;
		rotation[2] /*[0][2]*/ *= invInorm;

		rotation[3] /*[1][0]*/ *= invJnorm;
		rotation[4] /*[1][1]*/ *= invJnorm;
		rotation[5] /*[1][2]*/ *= invJnorm;

		/* row2 = row0 x row1 (cross product) */
		rotation[6] /*->m[2][0]*/ = rotation[1] /*->m[0][1]*/ * rotation[5] /*->m[1][2]*/ -
			rotation[2] /*->m[0][2]*/ * rotation[4] /*->m[1][1]*/;

		rotation[7] /*->m[2][1]*/ = rotation[2] /*->m[0][2]*/ * rotation[3] /*->m[1][0]*/ -
			rotation[0] /*->m[0][0]*/ * rotation[5] /*->m[1][2]*/;

		rotation[8] /*->m[2][2]*/ = rotation[0] /*->m[0][0]*/ * rotation[4] /*->m[1][1]*/ -
			rotation[1] /*->m[0][1]*/ * rotation[3] /*->m[1][0]*/;

		scale = (inorm + jnorm) / 2.0f;//1/s
		inv_Z = scale * inv_focalLength / cos(theta_M0) * g_M0;//1/Z0

		count++;
		pObject->ite_num = count;
		converged = ((criteria.type & CV_TERMCRIT_EPS) && (diff < criteria.epsilon));
		converged |= ((criteria.type & CV_TERMCRIT_ITER) && (count == criteria.max_iter));
	}
	const float invScale = 1 / scale;
	//const float invScale = 1 / (inv_Z * focalLength);//对于鱼眼用Z0/f作为平移向量的放大系数更可靠，注意此处invScale并不是代表取倒数，只是变量名没有改
	translation[0] = (imagePoints[0].x-xCenter) * invScale;
	translation[1] = (imagePoints[0].y-yCenter) * invScale;
	translation[2] = 1 / inv_Z;

	return CV_NO_ERR;
}


int  icvReleasePOSITObject(POSITObject ** ppObject)
{
	cvFree(ppObject);
	return CV_NO_ERR;
}

void cvPOSIT(POSITObject * pObject, CvPoint2D32f * imagePoints, double focalLength, CvTermCriteria criteria, float xCenter, float yCenter, int lensType, float* rotation, float* translation)
{
	icvPOSIT(pObject, imagePoints, (float)focalLength, criteria, xCenter, yCenter,
		lensType, rotation, translation);
}
POSITObject* CreatePOSITObject(CvPoint3D32f * points, int numPoints)
{
	POSITObject *pObject = 0;
	icvCreatePOSITObject(points, numPoints, &pObject);
	return pObject;
}

void cvReleasePOSITObject(POSITObject ** ppObject)
{
	icvReleasePOSITObject(ppObject);
}
