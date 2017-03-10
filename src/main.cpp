#include "absolute_orientation.h"
#include<iostream>
using namespace std;
int main(int argc, char **argv)
{
if(argc != 3)
    {
        cout << "Usage:  ./Coordinate_Align path_to_dataset_of_keyframe  path_to_dataset_of_groundtruth " << endl;
        return 1;
    }
    Coordinate_Align align;
    align.input_datas(argv[1],argv[2]);
    return 0;
}
