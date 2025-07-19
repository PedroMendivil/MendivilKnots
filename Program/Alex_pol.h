//Structs

struct knot_type
{
  int invariant;
  char name[8];
};

struct crossing_info
{
  int hand;
  int i_o;
  int i_u;
  float t_o;
  float t_u;
  int i;
  int j;
  int k;
};

typedef struct crossing_info crossing_info;
typedef struct knot_type knot_type;

const knot_type known_knot[] = 
{
  [0].invariant=      905463, [0].name="3_1",
  [1].invariant=     2509099, [1].name="4_1",
  [2].invariant=     2545745, [2].name="5_1",
  [3].invariant=     4925488, [3].name="5_2",
  [4].invariant=     8132760, [4].name="6_1",
  [5].invariant=    12240588, [5].name="6_2",
  [6].invariant=    17066075, [6].name="6_3",
  [7].invariant=     5080627, [7].name="7_1",
  [8].invariant=    12160074, [8].name="7_2",
  [9].invariant=    17161433, [9].name="7_3",
  [10].invariant=    22609223, [10].name="7_4",
  [11].invariant=    29272665, [11].name="7_5",
  [12].invariant=    36411894, [12].name="7_6",
  [13].invariant=    44444654, [13].name="7_7",
  [14].invariant=    16970983, [14].name="8_1",
  [15].invariant=    29649143, [15].name="8_2",
  [16].invariant=    29023769, [16].name="8_3",
  [17].invariant=    36551120, [17].name="8_4",
  [18].invariant=    45062872, [18].name="8_5",
  [19].invariant=    53487839, [19].name="8_6",
  [20].invariant=    53996320, [20].name="8_7",
  [21].invariant=    63138814, [21].name="8_8",
  [22].invariant=    63691163, [22].name="8_9",
  [23].invariant=    74235537, [23].name="8_10",
  [24].invariant=    73639120, [24].name="8_11",
  [25].invariant=    84681481, [25].name="8_12",
  [26].invariant=    84893732, [26].name="8_13",
  [27].invariant=    97004963, [27].name="8_14",
  [28].invariant=   110104951, [28].name="8_15",
  [29].invariant=   124614840, [29].name="8_16",
  [30].invariant=   139135171, [30].name="8_17",
  [31].invariant=   205711733, [31].name="8_18",
  [32].invariant=      972667, [32].name="8_19",
  [33].invariant=     8198629, [33].name="8_20",
  [34].invariant=    22718960, [34].name="8_21",
  [35].invariant=     8602498, [35].name="9_1",
  [36].invariant=    22609223, [36].name="9_2",
  [37].invariant=    37253791, [37].name="9_3",
  [38].invariant=    44752528, [38].name="9_4",
  [39].invariant=    53151207, [39].name="9_5",
  [40].invariant=    74834358, [40].name="9_6",
  [41].invariant=    85106248, [41].name="9_7",
  [42].invariant=    97004963, [42].name="9_8",
  [43].invariant=    98376073, [43].name="9_9",
  [44].invariant=   110346957, [44].name="9_10",
  [45].invariant=   111076775, [45].name="9_11",
  [46].invariant=   123585369, [46].name="9_12",
  [47].invariant=   138589445, [47].name="9_13",
  [48].invariant=   138047253, [48].name="9_14",
  [49].invariant=   153380336, [49].name="9_15",
  [50].invariant=   155390394, [50].name="9_16",
  [51].invariant=   154813639, [51].name="9_17",
  [52].invariant=   170046495, [52].name="9_18",
  [53].invariant=   169445857, [53].name="9_19",
  [54].invariant=   170952185, [54].name="9_20",
  [55].invariant=   186389865, [55].name="9_21",
  [56].invariant=   187969550, [56].name="9_22",
  [57].invariant=   204718107, [57].name="9_23",
  [58].invariant=   205711733, [58].name="9_24",
  [59].invariant=   222958011, [59].name="9_25",
  [60].invariant=   224340024, [60].name="9_26",
  [61].invariant=   243685842, [61].name="9_27",
  [62].invariant=   263925059, [62].name="9_28",
  [63].invariant=   263925059, [63].name="9_29",
  [64].invariant=   284874514, [64].name="9_30",
  [65].invariant=   306724656, [65].name="9_31",
  [66].invariant=   353171871, [66].name="9_32",
  [67].invariant=   377343177, [67].name="9_33",
  [68].invariant=   482280901, [68].name="9_34",
  [69].invariant=    73244041, [69].name="9_35",
  [70].invariant=   139407199, [70].name="9_36",
  [71].invariant=   204059022, [71].name="9_37",
  [72].invariant=   328437926, [72].name="9_38",
  [73].invariant=   305108326, [73].name="9_39",
  [74].invariant=   570041235, [74].name="9_40",
  [75].invariant=   242245383, [75].name="9_41",
  [76].invariant=     4976778, [76].name="9_42",
  [77].invariant=    17449975, [77].name="9_43",
  [78].invariant=    29148084, [78].name="9_44",
  [79].invariant=    53319390, [79].name="9_45",
  [80].invariant=     8132760, [80].name="9_46",
  [81].invariant=    74434274, [81].name="9_47",
  [82].invariant=    73441448, [82].name="9_48",
  [83].invariant=    63322107, [83].name="9_49",
  [84].invariant=    29023769, [84].name="10_1",
  [85].invariant=    55369240, [85].name="10_2",
  [86].invariant=    62773025, [86].name="10_3",
  [87].invariant=    73837058, [87].name="10_4",
  [88].invariant=   113042166, [88].name="10_5",
  [89].invariant=   140227359, [89].name="10_6",
  [90].invariant=   186704696, [90].name="10_7",
  [91].invariant=    86390821, [91].name="10_8",
  [92].invariant=   157132370, [92].name="10_9",
  [93].invariant=   204388432, [93].name="10_10",
  [94].invariant=   187019793, [94].name="10_11",
  [95].invariant=   225380122, [95].name="10_12",
  [96].invariant=   282929039, [96].name="10_13",
  [97].invariant=   331376270, [97].name="10_14",
  [98].invariant=   188921713, [98].name="10_15",
  [99].invariant=   223302331, [99].name="10_16",
  [100].invariant=   173388341, [100].name="10_17",
  [101].invariant=   305511092, [101].name="10_18",
  [102].invariant=   265428501, [102].name="10_19",
  [103].invariant=   123841753, [103].name="10_20",
  [104].invariant=   207039303, [104].name="10_21",
  [105].invariant=   244769809, [105].name="10_22",
  [106].invariant=   354476575, [106].name="10_23",
  [107].invariant=   305511092, [107].name="10_24",
  [108].invariant=   430119188, [108].name="10_25",
  [109].invariant=   378691749, [109].name="10_26",
  [110].invariant=   512679691, [110].name="10_27",
  [111].invariant=   283705017, [111].name="10_28",
  [112].invariant=   402892359, [112].name="10_29",
  [113].invariant=   452933450, [113].name="10_30",
  [114].invariant=   328020316, [114].name="10_31",
  [115].invariant=   484312489, [115].name="10_32",
  [116].invariant=   426294598, [116].name="10_33",
  [117].invariant=   138318216, [117].name="10_34",
  [118].invariant=   241886750, [118].name="10_35",
  [119].invariant=   262425887, [119].name="10_36",
  [120].invariant=   283705017, [120].name="10_37",
  [121].invariant=   351437316, [121].name="10_38",
  [122].invariant=   379140448, [122].name="10_39",
  [123].invariant=   571698481, [123].name="10_40",
  [124].invariant=   511110381, [124].name="10_41",
  [125].invariant=   664390717, [125].name="10_42",
  [126].invariant=   540107174, [126].name="10_43",
  [127].invariant=   632186650, [127].name="10_44",
  [128].invariant=   801532507, [128].name="10_45",
  [129].invariant=    99763102, [129].name="10_46",
  [130].invariant=   173691998, [130].name="10_47",
  [131].invariant=   246954867, [131].name="10_48",
  [132].invariant=   355783685, [132].name="10_49",
  [133].invariant=   286436402, [133].name="10_50",
  [134].invariant=   456382897, [134].name="10_51",
  [135].invariant=   354476575, [135].name="10_52",
  [136].invariant=   538496402, [136].name="10_53",
  [137].invariant=   225380122, [137].name="10_54",
  [138].invariant=   375997010, [138].name="10_55",
  [139].invariant=   430119188, [139].name="10_56",
  [140].invariant=   633931832, [140].name="10_57",
  [141].invariant=   425818806, [141].name="10_58",
  [142].invariant=   570041235, [142].name="10_59",
  [143].invariant=   731354331, [143].name="10_60",
  [144].invariant=   111565395, [144].name="10_61",
  [145].invariant=   208716151, [145].name="10_62",
  [146].invariant=   328437926, [146].name="10_63",
  [147].invariant=   267326687, [147].name="10_64",
  [148].invariant=   403822455, [148].name="10_65",
  [149].invariant=   573358132, [149].name="10_66",
  [150].invariant=   400578102, [150].name="10_67",
  [151].invariant=   328020316, [151].name="10_68",
  [152].invariant=   766121168, [152].name="10_69",
  [153].invariant=   455394089, [153].name="10_70",
  [154].invariant=   600641665, [154].name="10_71",
  [155].invariant=   542256987, [155].name="10_72",
  [156].invariant=   697546628, [156].name="10_73",
  [157].invariant=   400578102, [157].name="10_74",
  [158].invariant=   664390717, [158].name="10_75",
  [159].invariant=   330956795, [159].name="10_76",
  [160].invariant=   403822455, [160].name="10_77",
  [161].invariant=   482787246, [161].name="10_78",
  [162].invariant=   381408425, [162].name="10_79",
  [163].invariant=   514251407, [163].name="10_80",
  [164].invariant=   731977836, [164].name="10_81",
  [165].invariant=   408027545, [165].name="10_82",
  [166].invariant=   699989471, [166].name="10_83",
  [167].invariant=   768681169, [167].name="10_84",
  [168].invariant=   334764675, [168].name="10_85",
  [169].invariant=   733855620, [169].name="10_86",
  [170].invariant=   666774847, [170].name="10_87",
  [171].invariant=  1032094805, [171].name="10_88",
  [172].invariant=   991857978, [172].name="10_89",
  [173].invariant=   602342779, [173].name="10_90",
  [174].invariant=   546589063, [174].name="10_91",
  [175].invariant=   804804748, [175].name="10_92",
  [176].invariant=   456875464, [176].name="10_93",
  [177].invariant=   517416399, [177].name="10_94",
  [178].invariant=   840587430, [178].name="10_95",
  [179].invariant=   874925245, [179].name="10_96",
  [180].invariant=   763565437, [180].name="10_97",
  [181].invariant=   666774847, [181].name="10_98",
  [182].invariant=   672175175, [182].name="10_99",
  [183].invariant=   434458686, [183].name="10_100",
  [184].invariant=   730102457, [184].name="10_101",
  [185].invariant=   541720353, [185].name="10_102",
  [186].invariant=   571698481, [186].name="10_103",
  [187].invariant=   607476095, [187].name="10_104",
  [188].invariant=   838577642, [188].name="10_105",
  [189].invariant=   576699794, [189].name="10_106",
  [190].invariant=   875607196, [190].name="10_107",
  [191].invariant=   404285798, [191].name="10_108",
  [192].invariant=   739520554, [192].name="10_109",
  [193].invariant=   698155554, [193].name="10_110",
  [194].invariant=   602908636, [194].name="10_111",
  [195].invariant=   776410202, [195].name="10_112",
  [196].invariant=  1249966536, [196].name="10_113",
  [197].invariant=   878343868, [197].name="10_114",
  [198].invariant=  1202239772, [198].name="10_115",
  [199].invariant=   924843234, [199].name="10_116",
  [200].invariant=  1076349691, [200].name="10_117",
  [201].invariant=   963710742, [201].name="10_118",
  [202].invariant=  1035065795, [202].name="10_119",
  [203].invariant=  1113614483, [203].name="10_120",
  [204].invariant=  1341218615, [204].name="10_121",
  [205].invariant=  1119019553, [205].name="10_122",
  [206].invariant=  1498319904, [206].name="10_123",
  [207].invariant=      148922, [207].name="10_124",
  [208].invariant=    12484463, [208].name="10_125",
  [209].invariant=    36831640, [209].name="10_126",
  [210].invariant=    85747332, [210].name="10_127",
  [211].invariant=    12730743, [211].name="10_128",
  [212].invariant=    63138814, [212].name="10_129",
  [213].invariant=    29272665, [213].name="10_130",
  [214].invariant=    97004963, [214].name="10_131",
  [215].invariant=     2545745, [215].name="10_132",
  [216].invariant=    36411894, [216].name="10_133",
  [217].invariant=    54507206, [217].name="10_134",
  [218].invariant=   138318216, [218].name="10_135",
  [219].invariant=    22718960, [219].name="10_136",
  [220].invariant=    62955787, [220].name="10_137",
  [221].invariant=   124872290, [221].name="10_138",
  [222].invariant=      777715, [222].name="10_139",
  [223].invariant=     8198629, [223].name="10_140",
  [224].invariant=    45062872, [224].name="10_141",
  [225].invariant=    23384986, [225].name="10_142",
  [226].invariant=    74235537, [226].name="10_143",
  [227].invariant=   153665944, [227].name="10_144",
  [228].invariant=      883662, [228].name="10_145",
  [229].invariant=   109863211, [229].name="10_146",
  [230].invariant=    73639120, [230].name="10_147",
  [231].invariant=    97689315, [231].name="10_148",
  [232].invariant=   170952185, [232].name="10_149",
  [233].invariant=    85747332, [233].name="10_150",
  [234].invariant=   187653654, [234].name="10_151",
  [235].invariant=    11757893, [235].name="10_152",
  [236].invariant=       79269, [236].name="10_153",
  [237].invariant=    16780738, [237].name="10_154",
  [238].invariant=    63691163, [238].name="10_155",
  [239].invariant=   124614840, [239].name="10_156",
  [240].invariant=   244045806, [240].name="10_157",
  [241].invariant=   205381258, [241].name="10_158",
  [242].invariant=   154526966, [242].name="10_159",
  [243].invariant=    45217741, [243].name="10_160",
  [244].invariant=     2436278, [244].name="10_161",
  [245].invariant=   123841753, [245].name="10_162",
  [246].invariant=   263925059, [246].name="10_163",
  [247].invariant=   204388432, [247].name="10_164",
  [248].invariant=   153380336, [248].name="10_165"
};

//Functions

double calc_det(int N, double **M)
{
  int *P=(int*)malloc((N+1)*sizeof(int));
  for( int i=0; i<(N+1); i++)
  {
    P[i]=i;
  }
  for( int j=0; j<N; j++)
  {
    int i_max;
    double M_max=0.0;
    for( int i=j; i<N; i++)
    {
      if( fabs(M[i][j])>M_max)
      { 
        i_max=i;
        M_max=fabs(M[i][j]);
      }
    }
    if( M_max<__FLT_MIN__){ return 0.0;}
    if( i_max!=j)
    {
      int P_aux;
      double *M_aux;
      P_aux=P[j];
      P[j]=P[i_max];
      P[i_max]=P_aux;
      M_aux=M[j];
      M[j]=M[i_max];
      M[i_max]=M_aux;
      P[N]++;
    }
    for( int i=j+1; i<N; i++)
    {
      M[i][j]/=M[j][j];
      for( int k=j+1; k<N; k++)
      {
        M[i][k]-=M[i][j]*M[j][k];
      }
    }
  }
  double det=1.0;
  for( int i=0; i<N; i++)
  {
    det*=M[i][i];
  }
  if( (P[N]-N)%2!=0)
  {
    det*=(-1.0);
  }
  return det;
}

void calc_invariant( int N, float *r, char *knot_name)
{
  //close chain
  float r_cm[3];
  r_cm[0]=r_cm[1]=r_cm[2]=0.0;
  for( int i_p=0; i_p<N; i_p++)
  {
    r_cm[0] += r[3*i_p+0];
    r_cm[1] += r[3*i_p+1];
    r_cm[2] += r[3*i_p+2];
  }
  r_cm[0] /= N;
  r_cm[1] /= N;
  r_cm[2] /= N;
  float dist_max=0.0;
  for( int i_p=0; i_p<N; i_p++)
  {
    float dist=0.0;
    for( int i_c=0; i_c<3; i_c++)
    {
      dist += (r[3*i_p+i_c]-r_cm[i_c])*(r[3*i_p+i_c]-r_cm[i_c]);
    }
    dist = sqrt(dist);
    if( dist>dist_max)
    {
      dist_max = dist;
    }
  }
  float *r_l=(float*)malloc(3*(N+4)*sizeof(float));
  for( int i_p=0; i_p<N; i_p++)
  {
    for( int i_c=0; i_c<3; i_c++)
    {
      r_l[3*i_p+i_c]=r[3*i_p+i_c];
    }
  }
  float norm_aux[3]={0.0,0.0,0.0};
  for( int i_c=0; i_c<3; i_c++)
  {
    r_l[3*(N)+i_c]=r[3*(N-1)+i_c]-r_cm[i_c];
    r_l[3*(N+2)+i_c]=r[3*(0)+i_c]-r_cm[i_c];
    norm_aux[0] += r_l[3*(N+0)+i_c]*r_l[3*(N+0)+i_c];
    norm_aux[2] += r_l[3*(N+2)+i_c]*r_l[3*(N+2)+i_c];
  }
  norm_aux[0]=sqrt(norm_aux[0]);
  norm_aux[2]=sqrt(norm_aux[2]);
  for( int i_c=0; i_c<3; i_c++)
  {
    r_l[3*(N+0)+i_c] = 1.5*dist_max*r_l[3*(N+0)+i_c]/norm_aux[0];
    r_l[3*(N+2)+i_c] = 1.5*dist_max*r_l[3*(N+2)+i_c]/norm_aux[2];
    r_l[3*(N+1)+i_c] = r_l[3*(N+0)+i_c]+r_l[3*(N+2)+i_c];
    norm_aux[1] += r_l[3*(N+1)+i_c]*r_l[3*(N+1)+i_c];
  }
  norm_aux[1]=sqrt(norm_aux[1]);
  for( int i_c=0; i_c<3; i_c++)
  {
    r_l[3*(N+1)+i_c] = 1.5*dist_max*r_l[3*(N+1)+i_c]/norm_aux[1];
    r_l[3*(N+0)+i_c] += r_cm[i_c];
    r_l[3*(N+1)+i_c] += r_cm[i_c];
    r_l[3*(N+2)+i_c] += r_cm[i_c];
    r_l[3*(N+3)+i_c] = r[3*(0)+i_c];
  }
  //find crossings
  int n_crossings=0;
  crossing_info *crossing=NULL;
  for( int i_l=0; i_l<(N+3); i_l++)
  {
    for( int j_l=i_l+2; j_l<(N+3); j_l++)
    {
      float link_i[3], link_j[3], sep_ij[3];
      for( int i_c=0; i_c<3; i_c++)
      {
        link_i[i_c] = r_l[3*(i_l+1)+i_c]-r_l[3*i_l+i_c];
        link_j[i_c] = r_l[3*(j_l+1)+i_c]-r_l[3*j_l+i_c];
        sep_ij[i_c] = r_l[3*j_l+i_c]-r_l[3*i_l+i_c];
      }
      float cross_prod = link_i[0]*link_j[1]-link_i[1]*link_j[0];
      float t_i = (sep_ij[0]*link_j[1]-sep_ij[1]*link_j[0])/cross_prod;
      float t_j = (sep_ij[0]*link_i[1]-sep_ij[1]*link_i[0])/cross_prod;
      if( (0.0<t_i)&&(t_i<1.0)&&(0.0<t_j)&&(t_j<1.0))
      {
        n_crossings++;
        crossing=(crossing_info*)realloc(crossing,n_crossings*sizeof(crossing_info));
        if( (r_l[3*i_l+2]+t_i*link_i[2])>(r_l[3*j_l+2]+t_j*link_j[2]))
        {
          crossing[n_crossings-1].hand = (+1)*((cross_prod>0)-(cross_prod<0));
          crossing[n_crossings-1].i_o = i_l;
          crossing[n_crossings-1].i_u = j_l;
          crossing[n_crossings-1].t_o = t_i;
          crossing[n_crossings-1].t_u = t_j;
        }
        else
        {
          crossing[n_crossings-1].hand = (-1)*((cross_prod>0)-(cross_prod<0));
          crossing[n_crossings-1].i_o = j_l;
          crossing[n_crossings-1].i_u = i_l;
          crossing[n_crossings-1].t_o = t_j;
          crossing[n_crossings-1].t_u = t_i;
        }
      }
    }
  }
  free(r_l);
  //locate last underpass
  int idx_last_u;
  float t_last_u=-1.0;
  for( int i_l=N+2; i_l>0; i_l--)
  {
    for( int l=0; l<n_crossings; l++)
    {
      if( (i_l==crossing[l].i_u)&&(crossing[l].t_u>t_last_u))
      {
        t_last_u=crossing[l].t_u;
        idx_last_u=l;
      }
    }
    if( t_last_u>0.0)
    {
      break;
    }
  }
  //number arcs
  int arc_idx=0;
  int *idx_l=(int*)malloc(n_crossings*sizeof(int));
  float *t_l=(float*)malloc(n_crossings*sizeof(int));
  int *o_l=(int*)malloc(n_crossings*sizeof(int));
  for( int i_l=0; i_l<(N+3); i_l++)
  {
    int n_l=0;
    for( int l=0; l<n_crossings; l++)
    {
      if( (crossing[l].i_o==i_l))
      {
        int j_l_new=n_l;
        for( int j_l=0; j_l<n_l; j_l++)
        {
          if( (crossing[l].t_o<t_l[j_l]))
          {
            j_l_new = j_l; break;
          }
        }
        for( int j_l=n_l; j_l>j_l_new; j_l--)
        {
          idx_l[j_l] = idx_l[j_l-1];
          t_l[j_l] = t_l[j_l-1];
          o_l[j_l] = o_l[j_l-1];
        }
        idx_l[j_l_new] = l;
        t_l[j_l_new] = crossing[l].t_o;
        o_l[j_l_new] = 1;
        n_l++;
      }
      if( (crossing[l].i_u==i_l))
      {
        int j_l_new=n_l;
        for( int j_l=0; j_l<n_l; j_l++)
        {
          if( (crossing[l].t_u<t_l[j_l]))
          {
            j_l_new = j_l; break;
          }
        }
        for( int j_l=n_l; j_l>j_l_new; j_l--)
        {
          idx_l[j_l] = idx_l[j_l-1];
          t_l[j_l] = t_l[j_l-1];
          o_l[j_l] = o_l[j_l-1];
        }
        idx_l[j_l_new] = l;
        t_l[j_l_new] = crossing[l].t_u;
        o_l[j_l_new] = 0;
        n_l++;
      }
    }
    for( int j_l=0; j_l<n_l; j_l++)
    {
      if( o_l[j_l])
      {
        crossing[idx_l[j_l]].i=arc_idx;
      }
      else
      {
        crossing[idx_l[j_l]].j=arc_idx;
        if( idx_l[j_l]!=idx_last_u){ arc_idx++;}else{ arc_idx=0;}
        crossing[idx_l[j_l]].k=arc_idx;
      }
    }
  }
  free(idx_l);
  free(t_l);
  free(o_l);
  //make Alexander matrix
  int ***A=(int***)malloc(n_crossings*sizeof(int**));
  for( int i=0; i<n_crossings; i++)
  {
    A[i]=(int**)malloc(n_crossings*sizeof(int*));
    for( int j=0; j<n_crossings; j++)
    {
      A[i][j]=(int*)calloc(2,sizeof(int));
    }
  }
  for( int l=0; l<n_crossings; l++)
  {
    if( crossing[l].hand==1)
    {
      A[l][crossing[l].i][0] += (+1);
      A[l][crossing[l].i][1] += (-1);
      A[l][crossing[l].j][0] += (-1);
      A[l][crossing[l].k][1] += (+1);
    }
    else
    {
      A[l][crossing[l].i][0] += (+1);
      A[l][crossing[l].i][1] += (-1);
      A[l][crossing[l].j][1] += (+1);
      A[l][crossing[l].k][0] += (-1);
    }
  }
  //calculate knot invariant
  int n_M=n_crossings-1;
  double **M=(double**)malloc(n_M*sizeof(double*));
  for( int i=0; i<n_M; i++)
  {
    M[i]=(double*)malloc(n_M*sizeof(double));
  }
  double fp_k_inv=1.0;
  double t=-1.1;
  for( int i=0; i<n_M; i++)
  {
    for( int j=0; j<n_M; j++)
    {
      M[i][j]=A[i][j][0]+t*A[i][j][1];
    }
  }
  fp_k_inv*=calc_det(n_M,M);
  for( int i=0; i<n_M; i++)
  {
    for( int j=0; j<n_M; j++)
    {
      M[i][j]=A[i][j][0]+(1.0/t)*A[i][j][1];
    }
  }
  fp_k_inv*=calc_det(n_M,M);
  //return knot name
  int int_k_inv=(int)round(100000.0*fp_k_inv);
  if( int_k_inv==0)
  {
    snprintf(knot_name,sizeof(knot_name),"f");
    return;
  }
  if( int_k_inv==100000)
  {
    snprintf(knot_name,sizeof(knot_name),"0");
    return;
  }
  for( int i=0; i<249; i++)
  {
    if(int_k_inv==known_knot[i].invariant)
    {
      snprintf(knot_name,sizeof(knot_name),"%s",known_knot[i].name);
      return;
    }
  }
  snprintf(knot_name,sizeof(knot_name),"(%d)",int_k_inv);
}
