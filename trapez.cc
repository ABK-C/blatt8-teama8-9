#include <cmath>
#include <iostream>
#include <vector>

// Testfunktionen als Funktor
class Pol1 {
public:
  double operator()(double x) { return 3 * x + 2; }
};

class Pol2 {
public:
 double operator()(double x) { return -2*x*x+3*x+1; }
};

class Gauss {
public:
  double operator()(double x) { return 1 / (sqrt(M_PI * 2)) * exp(-x * x / 2); }
};

// berechnet Werte nach Trapezformel von I_0 bis I_N
template <class Functor> std::vector<double> trapez(Functor f, double a, double b, int N) {
  std::vector<double> I(N + 1); // Feld mit N+1 Eintraegen
  const double h = (b - a);
  I[0] = (h / 2) * (f(a) + f(b));
  for (int k = 1; k <= N; ++k) {
    int n = pow(2, k);
    double sum = 0;
    for(int i = 1; i < n; i += 2){
      sum += f(a+i*(h/n));
    }
    I[k] = I[k - 1]/(2) + h*sum/n; // setze k-ten Wert im Feld
  }
  return I;
}

// berechnet die Richardsonextrapolation aus I(k-1)  und I(k)
double richardson(double Iprev, double I) { return I + (I-Iprev)/3; }

// berechet Naeherungen ueber das Romberg-Verfahren
// I: Ergebnis von trapez()
std::vector<std::vector<double>> romberg(std::vector<double> I) {
  const int N = I.size() - 1;
  std::vector<std::vector<double>> R(N + 1);
  for (int k = 0; k <= N; ++k){
    R[k].push_back(I[k]);
  }
  for (int n = 1; n <= N; n++){
    for (int k = 0; k < N-(n-1) ; ++k) {
      double extr = R[k+1][n-1] + (R[k+1][n-1]-R[k][n-1])/(pow(2,2*n)-1);
    R[k].push_back(extr);
    }
  }
  return R;
}


void testeAufgabe1() {
  Pol1 f;
  std::vector<double> I_f = trapez(f, 0, 3, 3);
  for (double tf : I_f) {
    std::cout << "A1: f:" << tf << " == " << 19.5 << ":" << (tf == 19.5 ? "ja" : "nein") << std::endl;
  }
  Pol2 g;
  std::cout << "A1: g(1) = 2 ?" << (g(1) == 2 ? "ja" : "nein") << std::endl;
  std::vector<double> I_g = trapez(g, 0, 3, 3);
  std::cout << "A1: g0:" << I_g[0] << " == " << -10.5 << ":" << (I_g[0] == -10.5 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g1:" << I_g[1] << " == " << -3.75 << ":" << (I_g[1] == -3.75 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g2:" << I_g[2] << " == " << -2.0625 << ":" << (I_g[2] == -2.0625 ? "ja" : "nein") << std::endl;
  double rich = richardson(I_g[0], I_g[1]);
  std::cout << "A1: Richardson : " << rich << " : " << (rich == -1.5 ? "ja " : "nein") << std::endl;
}

void testeAufgabe2() {
  Pol1 f;
  std::vector<std::vector<double>> Rf = romberg(trapez(f, 0, 3, 3));
  bool alle_richtig = true;
  int entries = 0;
  for (auto row : Rf) {
    for (double val : row) {
      alle_richtig &= val == 19.5;
      ++entries;
    }
  }
  std::cout << "A2: alle Eintraege für f sind 1.5:" << (alle_richtig ? "ja" : "nein") << std::endl;
  std::cout << "A2: korrekte Zahl an Einträgen:" << (entries == 10 ? "ja" : "nein") << std::endl;
  Pol2 g;
  std::vector<std::vector<double>> Rg = romberg(trapez(g, 0, 3, 3));
  std::cout << "A2: R[1][1] und R[2][1] für g gleich -1.5: " << ((Rg[1][1] == -1.5) && (Rg[2][1] == -1.5) ? " ja " : " nein") << std::endl;
}


int main() {
  // Testfunktion:
  Pol1 f;
  std::cout << "f(0) = " << f(0) << '\n';
  // berechne Trapezformel fuer f
  std::vector<double> tf = trapez(f, 0., 3., 3);
  std::cout
      << "#############################################################\n";
  // Ausgabe:
  std::cout << "Trapez:\n";
  for (unsigned int i = 0; i < tf.size(); ++i) { // Schleife ueber Werte im Feld
    std::cout << "I_" << i << " = " << tf[i] << std::endl;
  }
  std::cout << "Romberg:\n";
  std::vector<std::vector<double> > R = romberg(tf);
  for(int k = 0, l = R.size() ;  k < l ; ++k) {
    for(int n = 0, m = R[k].size() ;  n < m ; ++n) {
      std::cout << R[k][n] << " ";
    }
    std::cout << std::endl;
  }

  Pol2 g1;
   std::vector<double> tg1 = trapez(g1, 0., 3., 3);
   std::cout << "Trapez:\n";
   for (unsigned int i = 0; i < tg1.size(); ++i) { // Schleife ueber Werte im Feld
    std::cout << "I_" << i << " = " << tg1[i] << std::endl;
   }

   Gauss h1;
   std::vector<double> th1 = trapez(h1, 0., 3., 3);
   std::cout << "Trapez:\n";
   for (unsigned int i = 0; i < th1.size(); ++i) { // Schleife ueber Werte im Feld
    std::cout << "I_" << i << " = " << th1[i] << std::endl;
   }

   std::cout<<richardson(tf[0], tf[1])<<"  "<<richardson(tf[1], tf[2])<<"  "<<richardson(tf[2], tf[3])<<std::endl;

   std::cout<<richardson(tg1[0], tg1[1])<<"  "<<richardson(tg1[1], tg1[2])<<"  "<<richardson(tg1[2], tg1[3])<<std::endl;

   std::cout<<richardson(th1[0], th1[1])<<"  "<<richardson(th1[1], th1[2])<<"  "<<richardson(th1[2], th1[3])<<std::endl;

  testeAufgabe1();
  testeAufgabe2();
  
}
