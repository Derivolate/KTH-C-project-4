#ifndef POINT_HPP
#define POINT_HPP
template<typename T>
class Point {
    private:
        T x;
        T y;
    public:
        Point(T xx=0.0, T yy=0.0) :x(xx), y(yy) { }

        bool operator==(const Point& Q) const{
            if(x==Q.X() && y==Q.Y())
                return true;
            return false;
        }

        Point<T> operator+(const Point& Q) const{
            return Point(x+Q.x,y+Q.y);
        }

        Point<T> operator-() const {
            return Point(-x,-y);
        }

        T  X() const{
            return x;
        }

        T  Y() const{
            return y;
        }

        void  zero() {
            x = y = 0.0;
        }
};
#endif