namespace{
  // Most directions are implicitly defined
  template<class T, class U> T GetValAs(const U& x){return x;}

#ifdef OSCLIB_STAN
  // But not this one (since you do not want to do it by accident)
  template<> double GetValAs<double>(const stan::math::var& x){return x.val();}
#endif
}
