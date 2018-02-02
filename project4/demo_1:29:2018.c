int logical_point_index_to_point_index(int logical_idx, int *dims)
{
    int idx = logical_index[1]*dims[0] + logical_index[0];
    return 2*idx;
}

void mycode()
{
    int idx = logical_point_index_to_point_index(logical_idx, dims);
    F[idx+0]; // x-component for logical_idx
    F[idx+1]; // y-component for logical_idx
}