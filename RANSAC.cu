__global__ void RANSAC_Fundamental(const Point2Df *src, const Point2Df *dst,int pts_num, const int *rand_list, float inlier_threshold, int iterations, int *inliers, float *fundamental)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
	//如果超出了迭代次数则返回
    if(idx >= iterations) return;   
    int rand_idx[4];    //随机列表索引
    Point2Df _src[4];   //源影像像点
    Point2Df _dst[4];   //匹配像点
    //每个线程计算一个F矩阵
    float *F = &fundamental[idx*9];
    //获取每个线程负责的随机列表索引
    rand_idx[0] = rand_list[idx*4+0];    rand_idx[1] = rand_list[idx*4+1];
    rand_idx[2] = rand_list[idx*4+2];    rand_idx[3] = rand_list[idx*4+3];
	//去除重复像对
    if(rand_idx[0] == rand_idx[1]) return;    if(rand_idx[0] == rand_idx[2]) return;
    if(rand_idx[0] == rand_idx[3]) return;    if(rand_idx[1] == rand_idx[2]) return;
    if(rand_idx[1] == rand_idx[3]) return;    if(rand_idx[2] == rand_idx[3]) return;
    //通过随机列表索引获取求解基础矩阵所需的4个像对
    for(int i=0; i < 4; i++) 
	{
        _src[i].x = src[rand_idx[i]].x;
        _src[i].y = src[rand_idx[i]].y;
        _dst[i].x = dst[rand_idx[i]].x;
        _dst[i].y = dst[rand_idx[i]].y;
    }
    //对基础矩阵归一化以提高噪声的鲁棒性
    //8点法计算基础矩阵
    int ret = GetFundamental(_src, _dst, F);
    //计算基础矩阵的内围点个数
    inliers[idx] = EvalFundamental(src, dst, pts_num, F, inlier_threshold);
}
