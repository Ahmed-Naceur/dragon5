static struct {
    int_32 nrecor;
    int_32 ninput;
    int_32 maxlvl;
    int_32 nstack;
    int_32 ixrlst;
    int_32 ioulst;
    int_32 idblst;
    int_32 nobjet;
    char cparin[13];
    char cdatin[73];
} header ;

static struct {
    int_32 ilines;
    int_32 ilevel;
    int_32 irecor;
    int_32 maskck[nmaskc];
    int_32 ipacki[nmaskc];
    char cparin[13];
    char myreco[73];
} record1 ;

static struct {
    int_32 indlin;
    int_32 idatin;
    float_32 adatin;
    double ddatin;
    int_32 idclin;
    int_32 idefin;
    int_32 iusein;
    char cparin[13];
    char cdatin[73];
} record2 ;
