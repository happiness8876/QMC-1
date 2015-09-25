#ifndef MYVECTOR_HPP
#define MYVECTOR_HPP
template <class T> class MyVector
{
  public:
    MyVector (int L1, int L2=1, int L3=1) : _L1(L1), _L2(L2), _L3(L3), _W2(L1*L2) { _element.resize (L1*L2*L3, 0); }
    const T& operator() (int i) const { return _element [i]; }
    T& operator() (int i) { return _element [i]; }
    const T& operator() (int i, int j) const { return _element [i+j*_L1]; }
    T& operator() (int i, int j) { return _element [i+j*_L1]; }
    const T& operator() (int i, int j, int k) const { return _element [i+j*_L1+k*_W2]; }
    T& operator() (int i, int j, int k) { return _element [i+j*_L1+k*_W2]; }
    void reset () { std::fill(_element.begin(), _element.end(), 0); }

  private:
    int _L1, _L2, _L3, _W2;
    std::vector<T> _element;
};
#endif
