#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <bitset>
using namespace std;
void print_menu();
class MultiplicationLong{

public:
    MultiplicationLong(){};
    long long direct_multiplying(long long num1, long long num2){
        return num1*num2;
    }
    long long fraction(long long num, long long divider){
        return num%divider;
    }
    bool is_complex_multiplication(long long a, long long b) {
        
        return a == 0 || b == 0 ? false : (a > 100 || b > 100);
    }
    long long shift(long long num){
        long long i = 0;
        while (num > 10){
            num/=10;
            i++;
        }
        return i;
    }
    long long find_max_pow(long long num1, long long num2) {
        long long pow_1 = shift(num1);
        long long pow_2 = shift(num2);
        long long max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    virtual long long calculate(long long num1, long long num2){return 0;};
    long long* multiplyPolynomials(long long (*arr_arr)[3]) {
        long long* result = new long long[6]();
        for (int i = 0; i <= 4; ++i) {
            result[i] = 0;
        }

        for (int i = 0; i <= 2; ++i) {
            for (int j = 0; j <= 2; ++j) {
                result[i + j] += is_complex_multiplication(arr_arr[0][i], arr_arr[1][j]) ? calculate(arr_arr[0][i], arr_arr[1][j]) : (arr_arr[1][j] *  arr_arr[0][i]);
            }
        }
        return result;
    }

};
class Caracuba : public MultiplicationLong {
public:
    Caracuba(){};
    long long calculate(long long num1, long long num2){
        cout << "=========" << endl;
        cout << num1 << " " << num2 << endl;
        long long sum = 0;
        long long arr [2] {num1, num2};
        long long arr1[3] {};
        long long arr2[3] {};
        long long arr_arr[2][3] {arr1[3], arr2[3]};
        long long max_pow = find_max_pow(num1, num2);
        long long x = pow(10, max_pow);
        find_common(arr, arr_arr, max_pow);
        long long* ur = multiplyPolynomials(arr_arr);
        for (int i=0; i <= 3; ++i) {
            sum += ur[i] * pow(x, i);
        }
        delete[] ur;
        return sum;
    };
private:
    void find_common(long long *arr, long long (*arr_arr)[3], long long max_pow) {
        for (long long i =0; i<2; ++i) {
            long long right_side = pow(10, max_pow);
            long long num1_fraction = fraction(arr[i], right_side);
            long long left_side = arr[i]/right_side;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = left_side;
            arr_arr[i][2] = 0;
            cout << arr_arr[i][0] << " " << arr_arr[i][1] << " " << arr_arr[i][2] << endl;
        }
    }
    long long shift(long long num){
        long long i = 0;
        while (num > 100){
            num/=10;
            i++;
        }
        return i;
    }
    
    
};
class ToomCook : public MultiplicationLong{
public:
    ToomCook(){};
    long long calculate(long long num1, long long num2){
        long long arr[2] {num1, num2};
        long long arr1[3] {};
        long long arr2[3] {};
        long long arr_arr[2][3] {arr1[3], arr2[3]};
        long long log_b_1 = log(num1)/log(10);
        long long log_b_2 = log(num2)/log(10);
        long long k = log_b_1/3 > log_b_2/3 ? log_b_1/3 + 1 : log_b_2/3 +1;
        for (long long i =0; i<2; ++i) {
            long long second_pow = pow(10, k);
            long long num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            long long second_side = arr[i]%second_pow;
            arr[i]/=second_pow;
            long long first_side = arr[i];
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
        }
        
        return evaluating_values(arr_arr, k);
    };
    long long evaluating_values(long long (*arr_arr)[3], long long max_pow){
        //evaluation
        long long _p0 = arr_arr[0][0] +  arr_arr[0][2];
        long long p0 = arr_arr[0][0];
        long long p1 = _p0 + arr_arr[0][1];
        long long p_1 = _p0 - arr_arr[0][1] ;
        long long p_2 = (p_1 +arr_arr[0][2])*2 -  arr_arr[0][0];
        long long p_inf = arr_arr[0][2];
        
        long long _q0 = arr_arr[1][0] +  arr_arr[1][2];
        long long q0 = arr_arr[1][0];
        long long q1 = _q0 + arr_arr[1][1];
        long long q_1 = _q0 - arr_arr[1][1] ;
        long long q_2 = (q_1 +arr_arr[1][2])*2 -  arr_arr[1][0];
        long long q_inf = arr_arr[1][2];
        
        //Pointwise multiplication
        long long r0 = p0 * q0;
        long long r1 = p1 * q1;
        long long r_1 = p_1 * q_1;
        long long r_2 = p_2 * q_2;
        long long r_inf = p_inf * q_inf;
        
        //Interpolation
        long long _r0 = r0;
        long long _r4 = r_inf;
        long long _r3 = (r_2 - r1)/3;
        long long _r1 = (r1 - r_1)/2;
        long long _r2 = (r_1 - r0);
        _r3 = (_r2 - _r3)/2 + 2*r_inf;
        _r2 = _r2 + _r1 - _r4;
        _r1 = _r1 - _r3;
        long long *r_pointer1[5] = {&_r0, &_r1, &_r2, &_r3, &_r4};
        long long sum = 0;
        long power_of_ten = pow(10, max_pow);
        for(int i = 0; i < 5; ++i){
            long long x_i = static_cast<long long>(pow(power_of_ten, i));
            sum += *r_pointer1[i] * x_i;
        }
        return sum;
    }
    long long find_max_pow(long long num1, long long num2) {
        long long pow_1 = shift(num1);
        long long pow_2 = shift(num2);
        long long max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        
        return max_pow;
        
    }
};
class Strassen : public MultiplicationLong {
public:
    Strassen(){};
        void find_common(long long *arr, long long (*arr_arr)[3], long long max_pow) {
        for (long long i =0; i<2; ++i) {
            long long second_pow = pow(10, max_pow);
            long long first_pow = pow(10, max_pow*2);
            long long num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            long long second_side = arr[i];
            long long first_side = arr[i] / first_pow;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
            cout << arr_arr[i][0] << " " << arr_arr[i][1] << " " << arr_arr[i][2] << endl;
        }
    }
    
    long long calculate(long long num1, long long num2){
        cout << "=========" << endl;
        cout << num1 << " " << num2 << endl;
        long long sum = 0;
        long long max_pow = find_max_pow(num1, num2);
        long long arr [2] {num1, num2};
        long long arr1[3] = {};
        long long arr2[3] = {};
        long long arr_arr[2][3] {arr1[3], arr2[3]};
        max_pow = !is_even(max_pow) ? ++max_pow : max_pow;
        long long x = pow(10, max_pow);
        find_common(arr, arr_arr, max_pow);
        long long* ur = multiplyPolynomials(arr_arr);
        for (int i=0; i <= 3; ++i) {
            sum += ur[i] * pow(x, i);
        }
        return sum;
    };
    long long find_max_pow(long long num1, long long num2) {
        long long pow_1 = min_pow_10(num1);
        long long pow_2 = min_pow_10(num2);
        long long max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    
    bool is_even(long long num){
        return num % 2 == 0;
    }
    long long min_pow_10(long long num){
        long long pow_10 = 0;
        long long result = 1;
        while(result <= num){
            pow_10++;
            result = pow(10, pow_10);
            num/=10;
        }
        return pow_10;
    }
};

class Cook: public MultiplicationLong {
public:
    Cook(){};
    bitset<4> calculate_bit(long long num1, long long num2){
        bitset<4> binaryRepresentation(num2);
        bitset<4> binaryRepresentation1(1<<4);
        bitset<4> inversed = binaryRepresentation^(binaryRepresentation1);
        return inversed;
    }
};

class Solovey_Strassen{
public:
    Solovey_Strassen(){};
    double input(){
        long long num = 0, k;
        while (num == 0) {
            cout << "input odd number: ";
            cin >> num;
            if(!is_even(num)){
                num = 0;
            }
            
        }
        cout << "input accuracy: ";
        cin >> k;
        return calculate_probability(num, k);
        
    }
    double calculate_probability(long long num, long long k){
        for(long long i = 0; i< k; ++i){
            long long a = rand() % (num-3) + 2;
            if(gcd(num, a) > 1){
                return 0;
            }
            else if(is_compared(pow(a, (num-1)/2), (a/num), num)){
                cout << pow(a, (num-1)/2) << " " << (a/num) << endl;
                return 0;
            }
        }
        return 1 - pow(2, -k);
    }
    bool is_compared(long long d, long b, long long a){
            return ((d-b)%a==0);
        }
    bool is_even(long long num){
        return num % 2 != 0;
    }
    long long gcd(long long a, long long b){
            long long t;
            while (b != 0){
                t = b;
                b = a%b;
                a = t;
            }
            return a;
        }
};
    

int main() {
//    print_menu();
    long long a = 13442132132132;
    long long b = 25323213214;
    Caracuba caracuba_method;
    cout << "CARACUBA: " << endl;
    cout << caracuba_method.calculate(a, b) << endl;
    cout << caracuba_method.direct_multiplying(a, b) << endl;
//
//    ToomCook toom_cook_method;
//    cout << "Toom Cook: " << endl;
//    cout << toom_cook_method.calculate(a, b) << endl;
//    cout << toom_cook_method.direct_multiplying(a, b) << endl;
    
    
    Strassen strassen_method;
    cout << "STRASSEN: " << endl;
    cout << strassen_method.calculate(a, b) << endl;
    cout << strassen_method.direct_multiplying(a, b) << endl;
    
//
////    Solovey_Strassen solovey_strassen_method;
////    cout << "SOLOVEY STRASSEN: " << endl;
////    cout << solovey_strassen_method.input() << endl;
//
//    Cook cook_method;
//    cout << "COOK: " << endl;
//    cout << cook_method.calculate(43133, 1263) << endl;
//    cout << cook_method.direct_multiplying(43133, 1263) << endl;
//    int type;
//    cin >> type;
//
//    switch (type) {
//        case 1:{
//            Caracuba generator;
//            generator.Caracuba();
//            break;
//        }
//        case 2:{
//            QuadraticCongruentialMethod generator;
//            generator.quadraticCongruentialMethod();
//            break;
//        }
//        case 3:
//        {
//            FibonachiNumbersMethod generator;
//            generator.fibonachiNumbersMethod();
//            break;
//        }
//        case 4:
//        {
//            InverseCongruentialMethod generator;
//            generator.inverseCongruentialMethod();
//            break;
//        }
//        case 5:
//        {
//            UnionMethod generator;
//            generator.unionMethod();
//            break;
//        }
//        case 6:
//        {
//            SigmaMethod generator;
//            generator.sigmaMethod();
//            break;
//        }
//        case 7:
//        {
//            PolarMethod generator;
//            generator.polarMehod();
//            break;
//        }
//        case 8:{
//            RelationMethod generator;
//            generator.relationMethod();
//            break;
//
//        }
//        case 9:{
//            LogarithmMethod generator;
//            generator.logarithmMethod();
//            break;
//
//        }
//        case 10:{
//            ArensMethod generator;
//            generator.arensMethod();
//            break;
//
//        }
//        default:{
//            cout << "Wrong input!!!\nOnly 1-10 requires";
//        }
//    }
//
//    return 1;
    
}


void print_menu() {
    cout << "Enter generator type: " << endl;
    cout << "\t1 - Linear Congruential Generator" << endl;
    cout << "\t2 - Quadratic Congruential Generator" << endl;
    cout << "\t3 - Fibonacci Numbers Generator" << endl;
    cout << "\t4 - Inverse Congruent SequenceGenerator" << endl;
    cout << "\t5 - Union Method Generator" << endl;
    cout << "\t6 - Sigma Method Generator" << endl;
    cout << "\t7 - Polar Coordinate Method Generator" << endl;
    cout << "\t8 - Relation Method Generator" << endl;
    cout << "\t9 - Logarithm Method Generator" << endl;
    cout << "\t10 - Arens Method Generator\n" << endl;
    
}

