__global__ void RANSAC_Fundamental(const Point2Df *src, const Point2Df *dst,int pts_num, const int *rand_list, float inlier_threshold, int iterations, int *inliers, float *fundamental)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//��������˵��������򷵻�
    if(idx >= iterations) return;   
    int rand_idx[4];    //����б�����
    Point2Df _src[4];   //ԴӰ�����
    Point2Df _dst[4];   //ƥ�����
    //ÿ���̼߳���һ��F����
    float *F = &fundamental[idx*9];
    //��ȡÿ���̸߳��������б�����
    rand_idx[0] = rand_list[idx*4+0];    rand_idx[1] = rand_list[idx*4+1];
    rand_idx[2] = rand_list[idx*4+2];    rand_idx[3] = rand_list[idx*4+3];
	//ȥ���ظ����
    if(rand_idx[0] == rand_idx[1]) return;    if(rand_idx[0] == rand_idx[2]) return;
    if(rand_idx[0] == rand_idx[3]) return;    if(rand_idx[1] == rand_idx[2]) return;
    if(rand_idx[1] == rand_idx[3]) return;    if(rand_idx[2] == rand_idx[3]) return;
    //ͨ������б�������ȡ���������������4�����
    for(int i=0; i < 4; i++) 
	{
        _src[i].x = src[rand_idx[i]].x;
        _src[i].y = src[rand_idx[i]].y;
        _dst[i].x = dst[rand_idx[i]].x;
        _dst[i].y = dst[rand_idx[i]].y;
    }
    //�Ի��������һ�������������³����
    //8�㷨�����������
    int ret = GetFundamental(_src, _dst, F);
    //��������������Χ�����
    inliers[idx] = EvalFundamental(src, dst, pts_num, F, inlier_threshold);
}
