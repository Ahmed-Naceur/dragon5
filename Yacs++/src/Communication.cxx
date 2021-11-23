#include "Communication.hxx"
using namespace std;

Communication::Communication(): compo_(0) {}

Communication::Communication(const Communication& comm) : compo_(comm.compo_) {}

Communication::~Communication() {}

int Communication::initialize(void * compo) {
    //Initialize the connection with YACS.
    char ret[64];
    compo_ = compo;
    cerr << "Communication: Initialize the connection with YACS  compo = " << compo_ << endl;
    int info = cp_cd(compo_, ret);
    cerr << "Communication: Initialize info = " << info << " ret= " << ret << endl;
    return info;
}

int Communication::terminate() {
    //Close the connection with YACS.
    cerr << "Communication: Close the connection with YACS  compo = " << compo_ << endl;
    int info =  cp_fin(compo_, CP_ARRET);
    cerr << "Communication: Terminate info = " << info << endl;
    return info;

}

// send datastream

int Communication::send(const int iteration, const string portName, const int& val ) {
    int info = cp_een(compo_, CP_ITERATION, 0., iteration, (char*) portName.c_str(),
                      1, (int*) &val);
    return info;
}

int Communication::send(const int iteration, const string portName, const int size, IntPtrConst& tab ) {
    const float fUNUSED=0;
    float temps(fUNUSED);
    int info = cp_een(compo_, CP_ITERATION, temps, iteration, (char*) portName.c_str(),
                      size, (int*) &(tab[0]));
    return info;
}

int Communication::send(const int iteration, const string portName, const float& val ) {
    int info = cp_ere(compo_, CP_ITERATION, 0., iteration, (char*) portName.c_str(),
                      1, (float*) &val);
    return info;
}

int Communication::send(const int iteration, const string portName, const int size, FloatPtrConst& tab ) {
    const float fUNUSED=0;
    float temps(fUNUSED);
    int info = cp_ere(compo_, CP_ITERATION, temps, iteration, (char*) portName.c_str(),
                      size, (float*) &(tab[0]));
    return info;
}

int Communication::send(const int iteration, const string portName, const double& val ) {
    int info = cp_edb(compo_, CP_ITERATION, double(0.), iteration, (char*) portName.c_str(),
                      1, (double*) &val);
    return info;
}

int Communication::send(const int iteration, const string portName, const int size, DoublePtrConst& tab ) {
    const float fUNUSED=0;
    double temps(fUNUSED);
    int info = cp_edb(compo_, CP_ITERATION, temps, iteration, (char*) portName.c_str(),
                      size, (double*) &(tab[0]));
    return info;
}

int Communication::send(const int iteration, const string portName, const bool& val ) {
    const float fUNUSED=0;
    float temps(fUNUSED);
    int* var;
    var = new int[1];
    if (val) {
      var[0] = 1;
    } else {
      var[0] = 0;
    }
    int info = cp_elo(compo_, CP_ITERATION, temps, iteration, (char*) portName.c_str(), 1, var);
    delete [] var;
    return info;
}

int Communication::send(const int iteration, const string portName, const int size, BoolPtrConst& tab ) {
    const float fUNUSED=0;
    float temps(fUNUSED);
    int* var;
    var = new int[size];
    for (int i=0; i < size; ++i) {
      if (tab[i]) {
        var[i] = 1;
      } else {
        var[i] = 0;
      }
    }
    int info = cp_elo(compo_, CP_ITERATION, temps, iteration, (char*) portName.c_str(), size, var);
    delete [] var;
    return info;
}

// receive datastream

int Communication::recv(int& iteration, const string portName, int& val ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    int* var;
    var = new int[1];
    int info = cp_len(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), 1, &nbValuesImported, var);
    val = var[0];
    delete [] var;
    return info;
}

int Communication::recv(int& iteration, const string portName, const int size, IntPtr& tab ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    int* var;
    var = new int[size];
    int info = cp_len(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), size, &nbValuesImported, var);
    for (int i=0; i < size; ++i) tab[i] = var[i];
    delete [] var;
    return info;    
}

int Communication::recv(int& iteration, const string portName, float& val ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    float* var;
    var = new float[1];
    int info = cp_lre(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), 1, &nbValuesImported, var);
    val = var[0];
    delete [] var;
    return info;
}

int Communication::recv(int& iteration, const string portName, const int size, FloatPtr& tab ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    float* var;
    var = new float[size];
    int info = cp_lre(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), size, &nbValuesImported, var);
    for (int i=0; i < size; ++i) tab[i] = var[i];
    delete [] var;
    return info;    
}

int Communication::recv(int& iteration, const string portName, double& val ) {
    double tempsi(0.) ;
    double tempsf(1.) ;
    int nbValuesImported;
    double* var;
    var = new double[1];
    int info = cp_ldb(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), 1, &nbValuesImported, var);
    val = var[0];
    delete [] var;
    return info;
}

int Communication::recv(int& iteration, const string portName, const int size, DoublePtr& tab ) {
    double tempsi(0.) ;
    double tempsf(1.) ;
    int nbValuesImported;
    double* var;
    var = new double[size];
    int info = cp_ldb(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), size, &nbValuesImported, var);
    for (int i=0; i < size; ++i) tab[i] = var[i];
    delete [] var;
    return info;    
}

int Communication::recv(int& iteration, const string portName, bool& val ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    int* var;
    var = new int[1];
    int info = cp_llo(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), 1, &nbValuesImported, (int*) var);
    val = (var[0] == 1);
    delete [] var;
    return info;
}

int Communication::recv(int& iteration, const string portName, const int size, BoolPtr& tab ) {
    float tempsi(0.) ;
    float tempsf(1.) ;
    int nbValuesImported;
    int* var;
    var = new int[size];
    int info = cp_llo(compo_, CP_ITERATION, &tempsi, &tempsf, &iteration, 
                      (char*) portName.c_str(), size, &nbValuesImported, (int*) var);
    for (int i=0; i < size; ++i) tab[i] = (var[i] == 1);
    delete [] var;
    return info;    
}
