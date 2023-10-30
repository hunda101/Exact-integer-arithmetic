#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <bitset>
#include <cassert>
#include <vector>

using namespace std;
void print_menu();
class BigInt {
private:
    vector<int> digits;
    
public:
    BigInt(){};
    BigInt(const std::string& numberStr) {
        for (char digit : numberStr) {
            if (isdigit(digit)) {
                digits.push_back(digit - '0');
            }
        }
        
        std::reverse(digits.begin(), digits.end());
    }
    BigInt(int number) {
            

            while (number > 0) {
                digits.push_back(number % 10);
                number /= 10;
            }

            if (digits.empty()) {
                digits.push_back(0);
            }
        }
    friend std::ostream& operator<<(std::ostream& os, const BigInt& bigInt) {

           if (bigInt.digits.empty()) {
               os << '0';
           } else {
               for (auto it = bigInt.digits.rbegin(); it != bigInt.digits.rend(); ++it) {
                   os << *it;
               }
           }

           return os;
       }
    BigInt big_pow(BigInt exponent) const {
        
        
        if (exponent == 0) {
            return 1;
        }
        
        BigInt result = *this;
        for (BigInt i = 1; i < exponent; i += 1) {
            result *= *this;
        }
        return result;
    }
    double log() const {

            double result = 0.0;
            for (size_t i = 0; i < digits.size(); ++i) {
                result = result * 10 + digits[i];
            }

            return std::log(result);
        }
    std::vector<int>::const_iterator begin() const {
            return digits.begin();
        }

        std::vector<int>::const_iterator end() const {
            return digits.end();
        }
    BigInt& operator+=(const BigInt& other) {
            
                int carry = 0;
                size_t max_size = std::max(digits.size(), other.digits.size());

                for (size_t i = 0; i < max_size || carry; ++i) {
                    if (i == digits.size()) {
                        digits.push_back(0);
                    }
                    int sum = digits[i] + carry + (i < other.digits.size() ? other.digits[i] : 0);
                    digits[i] = sum % 10;
                    carry = sum / 10;
                }
            

            return *this;
        }
    
    BigInt operator+(const BigInt& other) const {
        BigInt result = *this;
            result += other;
            return result;
        }
    BigInt operator+(int value) const {
        BigInt result = *this;
        result += BigInt(std::to_string(value));
        return result;
    }
    BigInt operator*(int value) const {
        BigInt result = *this;
        result *= BigInt(std::to_string(value));
        return result;
    }
    friend BigInt operator+(int value, const BigInt& bigInt) {
        return bigInt + value;
    }
    friend BigInt operator*(int value, const BigInt& bigInt) {
        return bigInt * value;
    }
    BigInt operator-(const BigInt& other) const {
            BigInt result = *this;
            result -= other;
            return result;
        }
    BigInt& operator-=(const BigInt& other) {
           
                if (*this < other) {
                    *this = other - *this;
                } else {
                    int borrow = 0;

                    for (size_t i = 0; i < digits.size(); ++i) {
                        int diff = digits[i] - borrow - (i < other.digits.size() ? other.digits[i] : 0);
                        if (diff < 0) {
                            diff += 10;
                            borrow = 1;
                        } else {
                            borrow = 0;
                        }
                        digits[i] = diff;
                    }

                    while (digits.size() > 1 && digits.back() == 0) {
                        digits.pop_back();
                    
                }
            }

            return *this;
        }
    
    BigInt operator/(const BigInt& other) const {
            if (other.isZero()) {
                throw std::runtime_error("Division by zero");
            }

            BigInt quotient;
            BigInt remainder = *this;

            while (remainder >= other) {
                int factor = 1;
                BigInt temp = other;

                while (remainder >= temp) {
                    remainder -= temp;
                    quotient += factor;
                    temp *= 10;
                    factor *= 10;
                }
            }

            return quotient;
        }
    BigInt& operator/=(const BigInt& other) {
            *this = *this / other;
            return *this;
        }
    BigInt operator*(const BigInt& other) const {
            BigInt result;
            result.digits.resize(digits.size() + other.digits.size(), 0);

            for (size_t i = 0; i < digits.size(); ++i) {
                int carry = 0;
                for (size_t j = 0; j < other.digits.size() || carry; ++j) {
                    long long cur = result.digits[i + j] + 1ll * digits[i] * (j < other.digits.size() ? other.digits[j] : 0) + carry;
                    result.digits[i + j] = static_cast<int>(cur % 10);
                    carry = static_cast<int>(cur / 10);
                }
            }

            while (result.digits.size() > 1 && result.digits.back() == 0) {
                result.digits.pop_back();
            }
            
            
            return result;
        }
    
    BigInt& operator*=(const BigInt& other) {
            *this = *this * other;
            return *this;
        }
    BigInt operator%(const BigInt& other) const {
        BigInt result = *this;
            result -= (result / other) * other;
            return result;
        }
    bool operator==(const BigInt& other) const {
            return digits == other.digits;
        }
    bool operator>(const BigInt& other) const {
        
            return digits.size() > other.digits.size() ||
                   (digits.size() == other.digits.size() && std::lexicographical_compare(digits.rbegin(), digits.rend(), other.digits.rbegin(), other.digits.rend()));
        
    }
    bool operator>=(const BigInt& other) const {
        return !(*this < other);
    }
    bool operator<(const BigInt& other) const {
    
            return digits.size() < other.digits.size() ||
                   (digits.size() == other.digits.size() && std::lexicographical_compare(digits.rbegin(), digits.rend(), other.digits.rbegin(), other.digits.rend()));
        
    }
    bool operator<=(const BigInt& other) const {
        return *this < other || *this == other;
    }
    private:
        bool isZero() const {
            return digits.size() == 1 && digits[0] == 0;
        }
    };
class MultiplicationLong{

public:
    MultiplicationLong(){};
    BigInt direct_multiplying(BigInt num1, BigInt num2){
        return num1*num2;
    }
    BigInt fraction(BigInt num, BigInt divider){
        return num%divider;
    }
    bool is_complex_multiplication(BigInt a, BigInt b) {
        
        return a ==  BigInt("0") || b ==  BigInt("0") ? false : (a > BigInt("10000") || b > BigInt("10000"));
    }
    BigInt shift(BigInt num){
        BigInt i = 0;
        BigInt ten = 10;
        while (num > 10){
            
            num= num / ten;
            i+=1;
        }
        
        return i;
    }
    BigInt find_max_pow(BigInt num1, BigInt num2) {
        
        BigInt pow_1 = shift(num1);
        BigInt pow_2 = shift(num2);
        
        BigInt max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    virtual BigInt calculate(BigInt num1, BigInt num2){return 0;};
    BigInt* multiplyPolynomials(BigInt (*arr_arr)[3]) {
        BigInt* result = new BigInt[6]();
        for (int i = 0; i <= 4; ++i) {
            result[i] = 0;
        }

        for (int i = 0; i <= 2; ++i) {
            for (int j = 0; j <= 2; ++j) {
                result[i + j] += is_complex_multiplication(arr_arr[1][j], arr_arr[0][i]) ? calculate(arr_arr[1][j], arr_arr[0][i]) : (arr_arr[1][j] *  arr_arr[0][i]);
            }
        }
        return result;
    }

};
class Caracuba : public MultiplicationLong {
public:
    Caracuba(){};
    BigInt calculate(BigInt num1, BigInt num2){
        BigInt sum = 0;
        BigInt arr [2] {num1, num2};
        BigInt arr1[3] {};
        BigInt arr2[3] {};
        
        BigInt arr_arr[2][3] {arr1[3], arr2[3]};
        
        BigInt max_pow = find_max_pow(num1, num2);
        
        BigInt base = 10;
        BigInt x = base.big_pow(max_pow);
        find_common(arr, arr_arr, max_pow);
        BigInt* ur = multiplyPolynomials(arr_arr);
        
        for (int i=0; i <= 3; ++i) {
            sum += ur[i] * x.big_pow(i);
        }
        delete[] ur;
        return sum;
    };
private:
    void find_common(BigInt *arr, BigInt (*arr_arr)[3], BigInt max_pow) {
        for (int i =0; i<2; ++i) {
            BigInt right_side = BigInt("10").big_pow(max_pow);
            BigInt num1_fraction = fraction(arr[i], right_side);
            BigInt left_side = arr[i]/right_side;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = left_side;
            arr_arr[i][2] = 0;
         }
    }
    BigInt shift(BigInt num){
        BigInt i = 0;
        while (num > 100){
            num/=10;
            i += 1;
        }
        return i;
    }
    
    
};
class ToomCook : public MultiplicationLong{
public:
    ToomCook(){};
    BigInt calculate(BigInt num1, BigInt num2){
        BigInt arr[2] {num1, num2};
        BigInt arr1[3] {};
        BigInt arr2[3] {};
        BigInt arr_arr[2][3] {arr1[3], arr2[3]};
        BigInt log_b_1 = num1.log()/log(10);
        BigInt log_b_2 = num2.log()/log(10);
        BigInt k = log_b_1/3 > log_b_2/3 ? log_b_1/3 + 1 : log_b_2/3 +1;
        for (int i =0; i<2; ++i) {
            BigInt second_pow = BigInt("10").big_pow(k);
            BigInt num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            BigInt second_side = arr[i]%second_pow;
            arr[i]/=second_pow;
            BigInt first_side = arr[i];
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
        }
        
        return evaluating_values(arr_arr, k);
    };
    BigInt evaluating_values(BigInt (*arr_arr)[3], BigInt max_pow){
        //evaluation
        BigInt _p0 = arr_arr[0][0] +  arr_arr[0][2];
        BigInt p0 = arr_arr[0][0];
        BigInt p1 = _p0 + arr_arr[0][1];
        BigInt p_1 = _p0 - arr_arr[0][1] ;
        BigInt p_2 = (p_1 +arr_arr[0][2])*2 -  arr_arr[0][0];
        BigInt p_inf = arr_arr[0][2];
        
        BigInt _q0 = arr_arr[1][0] +  arr_arr[1][2];
        BigInt q0 = arr_arr[1][0];
        BigInt q1 = _q0 + arr_arr[1][1];
        BigInt q_1 = _q0 - arr_arr[1][1] ;
        BigInt q_2 = (q_1 +arr_arr[1][2])*2 -  arr_arr[1][0];
        BigInt q_inf = arr_arr[1][2];
        
        //Pointwise multiplication
        BigInt r0 = p0 * q0;
        BigInt r1 = p1 * q1;
        BigInt r_1 = p_1 * q_1;
        BigInt r_2 = p_2 * q_2;
        BigInt r_inf = p_inf * q_inf;
        
        //Interpolation
        BigInt _r0 = r0;
        BigInt _r4 = r_inf;
        BigInt _r3 = (r_2 - r1)/3;
        BigInt _r1 = (r1 - r_1)/2;
        BigInt _r2 = (r_1 - r0);
        _r3 = (_r2 - _r3)/2 + 2*r_inf;
        _r2 = _r2 + _r1 - _r4;
        _r1 = _r1 - _r3;
        BigInt *r_pointer1[5] = {&_r0, &_r1, &_r2, &_r3, &_r4};
        BigInt sum = 0;
        BigInt power_of_ten = BigInt("10").big_pow(max_pow);
        cout << power_of_ten << endl;
        for(int i = 0; i < 5; ++i){
            BigInt x_i = power_of_ten.big_pow(i);
            cout << x_i << endl;
            sum += *r_pointer1[i] * x_i;
        }
        return sum;
    }
    BigInt find_max_pow(BigInt num1, BigInt num2) {
        BigInt pow_1 = shift(num1);
        BigInt pow_2 = shift(num2);
        BigInt max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        
        return max_pow;
        
    }
};
class Strassen : public MultiplicationLong {
public:
    Strassen(){};
        void find_common(BigInt *arr, BigInt (*arr_arr)[3], BigInt max_pow) {
        for (int i =0; i<2; ++i) {
            BigInt second_pow = BigInt("10").big_pow(max_pow);
            BigInt first_pow = BigInt("10").big_pow(max_pow*2);
            BigInt num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            BigInt second_side = arr[i];
            BigInt first_side = arr[i] / first_pow;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
        }
    }
    
    BigInt calculate(BigInt num1, BigInt num2){
        BigInt sum = 0;
        BigInt max_pow = find_max_pow(num1, num2);
        BigInt arr [2] {num1, num2};
        BigInt arr1[3] = {};
        BigInt arr2[3] = {};
        BigInt arr_arr[2][3] {arr1[3], arr2[3]};
        max_pow = !is_even(max_pow) ? max_pow+=1 : max_pow;
        BigInt x = BigInt("10").big_pow(max_pow);
        find_common(arr, arr_arr, max_pow);
        BigInt* ur = multiplyPolynomials(arr_arr);
        for (int i=0; i <= 3; ++i) {
            sum += ur[i] * x.big_pow(i);
        }
        return sum;
    };
    BigInt find_max_pow(BigInt num1, BigInt num2) {
        BigInt pow_1 = min_pow_10(num1);
        BigInt pow_2 = min_pow_10(num2);
        BigInt max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    
    bool is_even(BigInt num){
        return num % 2 == 0;
    }
    BigInt min_pow_10(BigInt num){
        BigInt pow_10 = 0;
        BigInt result = 1;
        while(result <= num){
            pow_10+=1;
            result =BigInt("10").big_pow(pow_10);
            num/=10;
        }
        return pow_10;
    }
};

class Cook: public MultiplicationLong {
public:
    Cook(){};
    double calculate_bit(BigInt num1, BigInt num2){
        
        double result = 1.0 / 9.0;
        int* binary = reinterpret_cast<int*>(&result);
            bitset<32> bits(*binary);
            cout  << bits << endl;

        return result;
    }
};

class PrimeTestBase{
public:
    PrimeTestBase(){};
    bool is_compared(long long d, long b, long long a){
            return ((d-b)%a==0);
        }
    bool odd(long long num){
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
    long long jacobi(long long a, long long b){
        long long g;
        assert(odd(b));
        if(a >=b) a%=b;
        if(a == 0) return 0;
        if(a == 1) return 1;
        if(a < 0){
            if((b-1)/2 % 2 == 0){
                return jacobi(-a, b);
                
            }
            else{
                return -jacobi(a, b);
            }
        }
        if(a % 2 == 0){
            if(((b*b-1)/8) % 2 ==0){
                return +jacobi(a/2, b);
            }
            else{
                return -jacobi(a/2, b);
            }
        }
        g = gcd(a, b);
        assert(odd(a));
        if(g == a){
            return 0;
        }
        else if(g != 1){
            return jacobi(g, b)* jacobi(a/g, b);
        }
        else if(((a-1)*(b-1)/4)%2 == 0){
            return + jacobi(b, a);
        }
        else{
            return - jacobi(b, a);
        }
    }
    bool is_power_two(long long num){
        for(int i = 1; i <= num; i*=2){
            if( i == num) return true;
        }
        return false;
    }
        
};
class Solovey_Strassen : public PrimeTestBase{
public:
    Solovey_Strassen(){};
    double input(){
        long long num = 0, k;
        while (num == 0) {
            cout << "input odd number: ";
            cin >> num;
            if(!odd(num)){
                num = 0;
            }
            
        }
        cout << "input accuracy: ";
        cin >> k;
        return calculate_probability(num, k);
        
    }
    double calculate_probability(long long p, long long k){
        for(long long i = 0; i< k; ++i){
            long long a = rand() % (p-3) + 2;
            if(gcd(a, p) > 1){
                return 0;
            }
            long long j = static_cast<long long>(pow(a, (p-1)/2)) % p;
            long long Jacob = jacobi(a, p);
            long long r = (j - Jacob);
            r = r < 0 ? r += p : r;
            if((r % p) != p){
                return 1;
            }
        }
        return 1 - pow(2, -k);
    }
    
};

class Luk_Lehmer : public PrimeTestBase{
public:
    Luk_Lehmer(){};
    double input(){
        long long p = 0;
        while (p == 0) {
            cout << "input odd number: ";
            cin >> p;
            if(!odd(p)){
                p = 0;
            }
            if (!(is_power_two(p+1))){
                p = 0;
            }
            
        }
        
        return calculate_probability(p);
        
    }
    double calculate_probability(long long p){
        long long S = 4;
        long long k = 1;
        long long M = static_cast<long long >(pow(2, p)) - 1;
        
        while(k != (p - 1)){
            S = ((S*S)- 2) % M;
            k += 1;
        }
        if(S == 0) return 1;
        return 0;
    }
    
};
class Rabin_Miller : public PrimeTestBase{
public:
    Rabin_Miller(){};
    double input(){
        long long p = 0, k = 0;
        while (p == 0) {
            cout << "input odd number: ";
            cin >> p;
            if(!odd(p)){
                p = 0;
            }
            
            
        }
        cout << "input accuracy: ";
        cin >> k;
        return calculate_probability(p, k);
        
    }
    double calculate_probability(long long p, long long k){
        
        int s = 0;
        long long m = p -1 ;
        while (m % 2 == 0) {
            m /= 2;
            s++;
        }
        for (int j = 0; j < k; ++j) {
            long long a = rand() % (p-3) + 2;
            long long x = static_cast<long>(pow(a, m))%p;
            long long y = (x*x)%p;
            for (int i = 0; i < s; ++i){
                if(y == 1 && x != 1 && x != (p-1)){
                    return 0;
                }
                x = y;
            }
            if( y != 1) return 0;
        }
        return 1;
        
    }
};
int main() {
//    print_menu();
    BigInt a("1051");
    BigInt b("2155");
    Caracuba caracuba_method;
    cout << "CARACUBA: " << endl;
    cout << caracuba_method.direct_multiplying(a, b) << endl;    cout << caracuba_method.calculate(a, b) << endl;
    ToomCook toom_cook_method;
    cout << "Toom Cook: " << endl;
    cout << toom_cook_method.calculate(a, b) << endl;
    cout << toom_cook_method.direct_multiplying(a, b) << endl;
    Strassen strassen_method;
    cout << "STRASSEN: " << endl;
    cout << strassen_method.calculate(a, b) << endl;
    cout << strassen_method.direct_multiplying(a, b) << endl;
//    Solovey_Strassen solovey_strassen_method;
//    cout << "SOLOVEY STRASSEN: " << endl;
//    cout << solovey_strassen_method.input() << endl;
    
    
//    Luk_Lehmer luk_lehmer;
//    cout << "luk_lehmer: " << endl;
//    cout << luk_lehmer.input() << endl;
    
    
    Rabin_Miller rabin_miller;
    cout << "rabin_miller: " << endl;
    cout << rabin_miller.input() << endl;
//
//    Cook cook_method;
//    cout << "COOK: " << endl;
//    cout << cook_method.calculate_bit(5, 9) << endl;
//    cout << cook_method.direct_multiplying(5, 9) << endl;
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

