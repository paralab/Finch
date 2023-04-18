//
// @author Milinda Fernando
// School of Computing, University of Utah

// NOTE: DO NOT CHANGE THIS FILE FOR ANY REASON.


// This header file contains all the rotation permutations + hilbert rotation table data hard corded to improve the performance of the hilbert curve.
// Note that: Rotations contains the concatenated strings of rot_perm and rot_index.
// Created by Milinda Fernando
// on 10/2/15.
//

#include "../include/hcurvedata.h"

char* rotations;
int* HILBERT_TABLE;

std::vector<unsigned int> RotationID_Stack;
unsigned int rotationStackPointer;


void _InitializeHcurve(int pDim) {

#ifdef HILBERT_ORDERING

if(pDim==2) {

    rotations = new char[_2D_ROTATIONS_SIZE];
    HILBERT_TABLE = new int[_2D_HILBERT_TABLE];



    strcpy(rotations + 0, "01230123\0");
    strcpy(rotations + 8, "03210321\0");
    strcpy(rotations + 16, "21032103\0");
    strcpy(rotations + 24, "23012301\0");



    /*HILBERT_TABLE[0]=1;
    HILBERT_TABLE[1]=0;
    HILBERT_TABLE[2]=0;
    HILBERT_TABLE[3]=2;
    HILBERT_TABLE[4]=0;
    HILBERT_TABLE[5]=1;
    HILBERT_TABLE[6]=1;
    HILBERT_TABLE[7]=3;
    HILBERT_TABLE[8]=3;
    HILBERT_TABLE[9]=2;
    HILBERT_TABLE[10]=2;
    HILBERT_TABLE[11]=0;
    HILBERT_TABLE[12]=2;
    HILBERT_TABLE[13]=3;
    HILBERT_TABLE[14]=3;
    HILBERT_TABLE[15]=1;*/


    HILBERT_TABLE[0] = 1;
    HILBERT_TABLE[1] = 0;
    HILBERT_TABLE[2] = 0;
    HILBERT_TABLE[3] = 2;
    HILBERT_TABLE[4] = 0;
    HILBERT_TABLE[5] = 3;
    HILBERT_TABLE[6] = 1;
    HILBERT_TABLE[7] = 1;
    HILBERT_TABLE[8] = 2;
    HILBERT_TABLE[9] = 2;
    HILBERT_TABLE[10] = 3;
    HILBERT_TABLE[11] = 0;
    HILBERT_TABLE[12] = 3;
    HILBERT_TABLE[13] = 1;
    HILBERT_TABLE[14] = 2;
    HILBERT_TABLE[15] = 3;

}

else if(pDim==3) {


    rotations = new char[_3D_ROTATIONS_SIZE];
    HILBERT_TABLE = new int[_3D_HILBERT_TABLE];

    RotationID_Stack.push_back(0);
    rotationStackPointer=0;


    /*
    //Table for  Canonical Hilbert index.

    strcpy(rotations + 0, "0123456701234567");
    strcpy(rotations + 16, "0743256107432561");
    strcpy(rotations + 32, "0167452301674523");
    strcpy(rotations + 48, "6107432521654703");
    strcpy(rotations + 64, "2543076147032165");
    strcpy(rotations + 80, "4523016745230167");
    strcpy(rotations + 96, "6125430761254307");
    strcpy(rotations + 112, "0761254303476521");
    strcpy(rotations + 128, "6701234523456701");
    strcpy(rotations + 144, "6745230167452301");
    strcpy(rotations + 160, "0347652107612543");
    strcpy(rotations + 176, "2103476521034765");
    strcpy(rotations + 192, "6547032147652103");
    strcpy(rotations + 208, "4567012345670123");
    strcpy(rotations + 224, "2165470361074325");
    strcpy(rotations + 240, "4307612525610743");
    strcpy(rotations + 256, "2561074343076125");
    strcpy(rotations + 272, "6521034743256107");
    strcpy(rotations + 288, "4703216525430761");
    strcpy(rotations + 304, "2345670167012345");
    strcpy(rotations + 320, "4325610765210347");
    strcpy(rotations + 336, "0321654703216547");
    strcpy(rotations + 352, "4765210365470321");
    strcpy(rotations + 368, "2301674523016745");

    HILBERT_TABLE[0] = 1;
    HILBERT_TABLE[1] = 2;
    HILBERT_TABLE[2] = 0;
    HILBERT_TABLE[3] = 3;
    HILBERT_TABLE[4] = 4;
    HILBERT_TABLE[5] = 0;
    HILBERT_TABLE[6] = 5;
    HILBERT_TABLE[7] = 6;
    HILBERT_TABLE[8] = 0;
    HILBERT_TABLE[9] = 9;
    HILBERT_TABLE[10] = 5;
    HILBERT_TABLE[11] = 8;
    HILBERT_TABLE[12] = 1;
    HILBERT_TABLE[13] = 1;
    HILBERT_TABLE[14] = 4;
    HILBERT_TABLE[15] = 7;
    HILBERT_TABLE[16] = 10;
    HILBERT_TABLE[17] = 0;
    HILBERT_TABLE[18] = 13;
    HILBERT_TABLE[19] = 14;
    HILBERT_TABLE[20] = 12;
    HILBERT_TABLE[21] = 2;
    HILBERT_TABLE[22] = 2;
    HILBERT_TABLE[23] = 11;
    HILBERT_TABLE[24] = 3;
    HILBERT_TABLE[25] = 6;
    HILBERT_TABLE[26] = 15;
    HILBERT_TABLE[27] = 3;
    HILBERT_TABLE[28] = 10;
    HILBERT_TABLE[29] = 11;
    HILBERT_TABLE[30] = 12;
    HILBERT_TABLE[31] = 14;
    HILBERT_TABLE[32] = 18;
    HILBERT_TABLE[33] = 12;
    HILBERT_TABLE[34] = 11;
    HILBERT_TABLE[35] = 17;
    HILBERT_TABLE[36] = 4;
    HILBERT_TABLE[37] = 16;
    HILBERT_TABLE[38] = 1;
    HILBERT_TABLE[39] = 4;
    HILBERT_TABLE[40] = 11;
    HILBERT_TABLE[41] = 5;
    HILBERT_TABLE[42] = 5;
    HILBERT_TABLE[43] = 12;
    HILBERT_TABLE[44] = 18;
    HILBERT_TABLE[45] = 13;
    HILBERT_TABLE[46] = 0;
    HILBERT_TABLE[47] = 17;
    HILBERT_TABLE[48] = 20;
    HILBERT_TABLE[49] = 3;
    HILBERT_TABLE[50] = 6;
    HILBERT_TABLE[51] = 6;
    HILBERT_TABLE[52] = 19;
    HILBERT_TABLE[53] = 2;
    HILBERT_TABLE[54] = 9;
    HILBERT_TABLE[55] = 0;
    HILBERT_TABLE[56] = 21;
    HILBERT_TABLE[57] = 18;
    HILBERT_TABLE[58] = 17;
    HILBERT_TABLE[59] = 22;
    HILBERT_TABLE[60] = 16;
    HILBERT_TABLE[61] = 7;
    HILBERT_TABLE[62] = 7;
    HILBERT_TABLE[63] = 1;
    HILBERT_TABLE[64] = 8;
    HILBERT_TABLE[65] = 22;
    HILBERT_TABLE[66] = 21;
    HILBERT_TABLE[67] = 8;
    HILBERT_TABLE[68] = 23;
    HILBERT_TABLE[69] = 18;
    HILBERT_TABLE[70] = 17;
    HILBERT_TABLE[71] = 9;
    HILBERT_TABLE[72] = 19;
    HILBERT_TABLE[73] = 1;
    HILBERT_TABLE[74] = 20;
    HILBERT_TABLE[75] = 9;
    HILBERT_TABLE[76] = 9;
    HILBERT_TABLE[77] = 7;
    HILBERT_TABLE[78] = 6;
    HILBERT_TABLE[79] = 8;
    HILBERT_TABLE[80] = 2;
    HILBERT_TABLE[81] = 19;
    HILBERT_TABLE[82] = 12;
    HILBERT_TABLE[83] = 21;
    HILBERT_TABLE[84] = 10;
    HILBERT_TABLE[85] = 10;
    HILBERT_TABLE[86] = 13;
    HILBERT_TABLE[87] = 23;
    HILBERT_TABLE[88] = 11;
    HILBERT_TABLE[89] = 14;
    HILBERT_TABLE[90] = 4;
    HILBERT_TABLE[91] = 6;
    HILBERT_TABLE[92] = 1;
    HILBERT_TABLE[93] = 3;
    HILBERT_TABLE[94] = 18;
    HILBERT_TABLE[95] = 11;
    HILBERT_TABLE[96] = 15;
    HILBERT_TABLE[97] = 4;
    HILBERT_TABLE[98] = 10;
    HILBERT_TABLE[99] = 12;
    HILBERT_TABLE[100] = 12;
    HILBERT_TABLE[101] = 17;
    HILBERT_TABLE[102] = 3;
    HILBERT_TABLE[103] = 16;
    HILBERT_TABLE[104] = 3;
    HILBERT_TABLE[105] = 13;
    HILBERT_TABLE[106] = 2;
    HILBERT_TABLE[107] = 16;
    HILBERT_TABLE[108] = 15;
    HILBERT_TABLE[109] = 5;
    HILBERT_TABLE[110] = 13;
    HILBERT_TABLE[111] = 4;
    HILBERT_TABLE[112] = 22;
    HILBERT_TABLE[113] = 11;
    HILBERT_TABLE[114] = 19;
    HILBERT_TABLE[115] = 2;
    HILBERT_TABLE[116] = 9;
    HILBERT_TABLE[117] = 0;
    HILBERT_TABLE[118] = 14;
    HILBERT_TABLE[119] = 14;
    HILBERT_TABLE[120] = 15;
    HILBERT_TABLE[121] = 15;
    HILBERT_TABLE[122] = 3;
    HILBERT_TABLE[123] = 20;
    HILBERT_TABLE[124] = 13;
    HILBERT_TABLE[125] = 23;
    HILBERT_TABLE[126] = 2;
    HILBERT_TABLE[127] = 19;
    HILBERT_TABLE[128] = 8;
    HILBERT_TABLE[129] = 5;
    HILBERT_TABLE[130] = 23;
    HILBERT_TABLE[131] = 13;
    HILBERT_TABLE[132] = 7;
    HILBERT_TABLE[133] = 4;
    HILBERT_TABLE[134] = 16;
    HILBERT_TABLE[135] = 16;
    HILBERT_TABLE[136] = 23;
    HILBERT_TABLE[137] = 13;
    HILBERT_TABLE[138] = 17;
    HILBERT_TABLE[139] = 17;
    HILBERT_TABLE[140] = 21;
    HILBERT_TABLE[141] = 12;
    HILBERT_TABLE[142] = 8;
    HILBERT_TABLE[143] = 5;
    HILBERT_TABLE[144] = 18;
    HILBERT_TABLE[145] = 18;
    HILBERT_TABLE[146] = 0;
    HILBERT_TABLE[147] = 9;
    HILBERT_TABLE[148] = 5;
    HILBERT_TABLE[149] = 8;
    HILBERT_TABLE[150] = 11;
    HILBERT_TABLE[151] = 22;
    HILBERT_TABLE[152] = 9;
    HILBERT_TABLE[153] = 10;
    HILBERT_TABLE[154] = 14;
    HILBERT_TABLE[155] = 23;
    HILBERT_TABLE[156] = 19;
    HILBERT_TABLE[157] = 21;
    HILBERT_TABLE[158] = 22;
    HILBERT_TABLE[159] = 19;
    HILBERT_TABLE[160] = 6;
    HILBERT_TABLE[161] = 20;
    HILBERT_TABLE[162] = 20;
    HILBERT_TABLE[163] = 15;
    HILBERT_TABLE[164] = 22;
    HILBERT_TABLE[165] = 10;
    HILBERT_TABLE[166] = 14;
    HILBERT_TABLE[167] = 21;
    HILBERT_TABLE[168] = 7;
    HILBERT_TABLE[169] = 15;
    HILBERT_TABLE[170] = 21;
    HILBERT_TABLE[171] = 10;
    HILBERT_TABLE[172] = 17;
    HILBERT_TABLE[173] = 21;
    HILBERT_TABLE[174] = 16;
    HILBERT_TABLE[175] = 20;
    HILBERT_TABLE[176] = 14;
    HILBERT_TABLE[177] = 22;
    HILBERT_TABLE[178] = 6;
    HILBERT_TABLE[179] = 7;
    HILBERT_TABLE[180] = 20;
    HILBERT_TABLE[181] = 1;
    HILBERT_TABLE[182] = 22;
    HILBERT_TABLE[183] = 18;
    HILBERT_TABLE[184] = 23;
    HILBERT_TABLE[185] = 20;
    HILBERT_TABLE[186] = 16;
    HILBERT_TABLE[187] = 19;
    HILBERT_TABLE[188] = 8;
    HILBERT_TABLE[189] = 15;
    HILBERT_TABLE[190] = 7;
    HILBERT_TABLE[191] = 23;

*/



    // Hilbert table based on the morton ordering.

    strcpy(rotations + 0, "0231576403127465\0");
    strcpy(rotations + 16, "0451376203741265\0");
    strcpy(rotations + 32, "0264573107163425\0");
    strcpy(rotations + 48, "6204513725163407\0");
    strcpy(rotations + 64, "3751046243705261\0");
    strcpy(rotations + 80, "5731026443527061\0");
    strcpy(rotations + 96, "6237510465127403\0");
    strcpy(rotations + 112, "0462375107341625\0");
    strcpy(rotations + 128, "6402315725341607\0");
    strcpy(rotations + 144, "6457310265741203\0");
    strcpy(rotations + 160, "0154673201763245\0");
    strcpy(rotations + 176, "3201546723105467\0");
    strcpy(rotations + 192, "6754013245763201\0");
    strcpy(rotations + 208, "5764023147563021\0");
    strcpy(rotations + 224, "3267540167105423\0");
    strcpy(rotations + 240, "5104623721563047\0");
    strcpy(rotations + 256, "3762045147305621\0");
    strcpy(rotations + 272, "6732015445327601\0");
    strcpy(rotations + 288, "5401326723541067\0");
    strcpy(rotations + 304, "3157640261705243\0");
    strcpy(rotations + 320, "5137620461527043\0");
    strcpy(rotations + 336, "0132675401327645\0");
    strcpy(rotations + 352, "5467320167541023\0");
    strcpy(rotations + 368, "3102645721305647\0");


    HILBERT_TABLE[0] = 1;
    HILBERT_TABLE[1] = 3;
    HILBERT_TABLE[2] = 2;
    HILBERT_TABLE[3] = 0;
    HILBERT_TABLE[4] = 6;
    HILBERT_TABLE[5] = 4;
    HILBERT_TABLE[6] = 5;
    HILBERT_TABLE[7] = 0;
    HILBERT_TABLE[8] = 0;
    HILBERT_TABLE[9] = 8;
    HILBERT_TABLE[10] = 9;
    HILBERT_TABLE[11] = 5;
    HILBERT_TABLE[12] = 7;
    HILBERT_TABLE[13] = 1;
    HILBERT_TABLE[14] = 4;
    HILBERT_TABLE[15] = 1;
    HILBERT_TABLE[16] = 10;
    HILBERT_TABLE[17] = 14;
    HILBERT_TABLE[18] = 0;
    HILBERT_TABLE[19] = 13;
    HILBERT_TABLE[20] = 11;
    HILBERT_TABLE[21] = 12;
    HILBERT_TABLE[22] = 2;
    HILBERT_TABLE[23] = 2;
    HILBERT_TABLE[24] = 3;
    HILBERT_TABLE[25] = 3;
    HILBERT_TABLE[26] = 6;
    HILBERT_TABLE[27] = 15;
    HILBERT_TABLE[28] = 14;
    HILBERT_TABLE[29] = 10;
    HILBERT_TABLE[30] = 12;
    HILBERT_TABLE[31] = 11;
    HILBERT_TABLE[32] = 18;
    HILBERT_TABLE[33] = 17;
    HILBERT_TABLE[34] = 12;
    HILBERT_TABLE[35] = 11;
    HILBERT_TABLE[36] = 4;
    HILBERT_TABLE[37] = 4;
    HILBERT_TABLE[38] = 1;
    HILBERT_TABLE[39] = 16;
    HILBERT_TABLE[40] = 11;
    HILBERT_TABLE[41] = 12;
    HILBERT_TABLE[42] = 5;
    HILBERT_TABLE[43] = 5;
    HILBERT_TABLE[44] = 17;
    HILBERT_TABLE[45] = 18;
    HILBERT_TABLE[46] = 0;
    HILBERT_TABLE[47] = 13;
    HILBERT_TABLE[48] = 20;
    HILBERT_TABLE[49] = 6;
    HILBERT_TABLE[50] = 3;
    HILBERT_TABLE[51] = 6;
    HILBERT_TABLE[52] = 0;
    HILBERT_TABLE[53] = 19;
    HILBERT_TABLE[54] = 9;
    HILBERT_TABLE[55] = 2;
    HILBERT_TABLE[56] = 21;
    HILBERT_TABLE[57] = 22;
    HILBERT_TABLE[58] = 18;
    HILBERT_TABLE[59] = 17;
    HILBERT_TABLE[60] = 1;
    HILBERT_TABLE[61] = 16;
    HILBERT_TABLE[62] = 7;
    HILBERT_TABLE[63] = 7;
    HILBERT_TABLE[64] = 8;
    HILBERT_TABLE[65] = 8;
    HILBERT_TABLE[66] = 22;
    HILBERT_TABLE[67] = 21;
    HILBERT_TABLE[68] = 9;
    HILBERT_TABLE[69] = 23;
    HILBERT_TABLE[70] = 17;
    HILBERT_TABLE[71] = 18;
    HILBERT_TABLE[72] = 19;
    HILBERT_TABLE[73] = 9;
    HILBERT_TABLE[74] = 1;
    HILBERT_TABLE[75] = 20;
    HILBERT_TABLE[76] = 8;
    HILBERT_TABLE[77] = 9;
    HILBERT_TABLE[78] = 6;
    HILBERT_TABLE[79] = 7;
    HILBERT_TABLE[80] = 2;
    HILBERT_TABLE[81] = 21;
    HILBERT_TABLE[82] = 19;
    HILBERT_TABLE[83] = 12;
    HILBERT_TABLE[84] = 23;
    HILBERT_TABLE[85] = 10;
    HILBERT_TABLE[86] = 13;
    HILBERT_TABLE[87] = 10;
    HILBERT_TABLE[88] = 11;
    HILBERT_TABLE[89] = 6;
    HILBERT_TABLE[90] = 14;
    HILBERT_TABLE[91] = 4;
    HILBERT_TABLE[92] = 11;
    HILBERT_TABLE[93] = 1;
    HILBERT_TABLE[94] = 18;
    HILBERT_TABLE[95] = 3;
    HILBERT_TABLE[96] = 15;
    HILBERT_TABLE[97] = 12;
    HILBERT_TABLE[98] = 4;
    HILBERT_TABLE[99] = 10;
    HILBERT_TABLE[100] = 16;
    HILBERT_TABLE[101] = 12;
    HILBERT_TABLE[102] = 3;
    HILBERT_TABLE[103] = 17;
    HILBERT_TABLE[104] = 3;
    HILBERT_TABLE[105] = 16;
    HILBERT_TABLE[106] = 13;
    HILBERT_TABLE[107] = 2;
    HILBERT_TABLE[108] = 4;
    HILBERT_TABLE[109] = 15;
    HILBERT_TABLE[110] = 13;
    HILBERT_TABLE[111] = 5;
    HILBERT_TABLE[112] = 22;
    HILBERT_TABLE[113] = 2;
    HILBERT_TABLE[114] = 11;
    HILBERT_TABLE[115] = 19;
    HILBERT_TABLE[116] = 14;
    HILBERT_TABLE[117] = 9;
    HILBERT_TABLE[118] = 14;
    HILBERT_TABLE[119] = 0;
    HILBERT_TABLE[120] = 15;
    HILBERT_TABLE[121] = 20;
    HILBERT_TABLE[122] = 15;
    HILBERT_TABLE[123] = 3;
    HILBERT_TABLE[124] = 19;
    HILBERT_TABLE[125] = 13;
    HILBERT_TABLE[126] = 2;
    HILBERT_TABLE[127] = 23;
    HILBERT_TABLE[128] = 8;
    HILBERT_TABLE[129] = 13;
    HILBERT_TABLE[130] = 5;
    HILBERT_TABLE[131] = 23;
    HILBERT_TABLE[132] = 16;
    HILBERT_TABLE[133] = 7;
    HILBERT_TABLE[134] = 16;
    HILBERT_TABLE[135] = 4;
    HILBERT_TABLE[136] = 23;
    HILBERT_TABLE[137] = 17;
    HILBERT_TABLE[138] = 13;
    HILBERT_TABLE[139] = 17;
    HILBERT_TABLE[140] = 5;
    HILBERT_TABLE[141] = 21;
    HILBERT_TABLE[142] = 8;
    HILBERT_TABLE[143] = 12;
    HILBERT_TABLE[144] = 18;
    HILBERT_TABLE[145] = 9;
    HILBERT_TABLE[146] = 18;
    HILBERT_TABLE[147] = 0;
    HILBERT_TABLE[148] = 22;
    HILBERT_TABLE[149] = 5;
    HILBERT_TABLE[150] = 11;
    HILBERT_TABLE[151] = 8;
    HILBERT_TABLE[152] = 9;
    HILBERT_TABLE[153] = 23;
    HILBERT_TABLE[154] = 10;
    HILBERT_TABLE[155] = 14;
    HILBERT_TABLE[156] = 19;
    HILBERT_TABLE[157] = 19;
    HILBERT_TABLE[158] = 22;
    HILBERT_TABLE[159] = 21;
    HILBERT_TABLE[160] = 6;
    HILBERT_TABLE[161] = 15;
    HILBERT_TABLE[162] = 20;
    HILBERT_TABLE[163] = 20;
    HILBERT_TABLE[164] = 21;
    HILBERT_TABLE[165] = 22;
    HILBERT_TABLE[166] = 14;
    HILBERT_TABLE[167] = 10;
    HILBERT_TABLE[168] = 7;
    HILBERT_TABLE[169] = 10;
    HILBERT_TABLE[170] = 15;
    HILBERT_TABLE[171] = 21;
    HILBERT_TABLE[172] = 20;
    HILBERT_TABLE[173] = 17;
    HILBERT_TABLE[174] = 16;
    HILBERT_TABLE[175] = 21;
    HILBERT_TABLE[176] = 14;
    HILBERT_TABLE[177] = 7;
    HILBERT_TABLE[178] = 22;
    HILBERT_TABLE[179] = 6;
    HILBERT_TABLE[180] = 18;
    HILBERT_TABLE[181] = 20;
    HILBERT_TABLE[182] = 22;
    HILBERT_TABLE[183] = 1;
    HILBERT_TABLE[184] = 23;
    HILBERT_TABLE[185] = 19;
    HILBERT_TABLE[186] = 20;
    HILBERT_TABLE[187] = 16;
    HILBERT_TABLE[188] = 23;
    HILBERT_TABLE[189] = 8;
    HILBERT_TABLE[190] = 7;
    HILBERT_TABLE[191] = 15;

    }


#else

const int num_orthant = (1u<<pDim);
rotations = new char[2*num_orthant];
HILBERT_TABLE = new int[num_orthant];
for (int ort = 0; ort < num_orthant; ort++)
{
  rotations[ort] = '0' + ort;
  rotations[ort + num_orthant] = '0' + ort;
  HILBERT_TABLE[ort] = 0;
}


/// if(pDim==2) {
/// 
///     rotations = new char[8];
///     HILBERT_TABLE = new char[4];
/// 
///     strcpy(rotations + 0, "0123401234\0");
/// 
///     HILBERT_TABLE[0] = 0;
///     HILBERT_TABLE[1] = 0;
///     HILBERT_TABLE[2] = 0;
///     HILBERT_TABLE[3] = 0;
/// }
/// 
/// else if(pDim==3) {
/// 
///     rotations = new char[16];
///     HILBERT_TABLE = new char[8];
/// 
///     memcpy(rotations + 0, "0123456701234567", 16);
/// 
///     HILBERT_TABLE[0] = 0;
///     HILBERT_TABLE[1] = 0;
///     HILBERT_TABLE[2] = 0;
///     HILBERT_TABLE[3] = 0;
///     HILBERT_TABLE[4] = 0;
///     HILBERT_TABLE[5] = 0;
///     HILBERT_TABLE[6] = 0;
///     HILBERT_TABLE[7] = 0;
/// 
///     }

#endif

}
