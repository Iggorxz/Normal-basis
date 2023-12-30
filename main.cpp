#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <chrono>
#include <unordered_map>

class GF2mElement {
private:
    std::vector<bool> coefficients;
    static const int m = 233;
    static const int p = 467;
    static std::vector<std::vector<int>> nonZeroColumns;
    static std::unordered_map<int, int> mod_pow_2_cache;

public:
    GF2mElement(const std::vector<bool>& coeffs) : coefficients(coeffs) {
        coefficients.resize(m, false);
    }

    GF2mElement(const std::string& bitString) {
        for (int i = bitString.length() - 1; i >= 0; --i) {
            coefficients.push_back(bitString[i] == '1');
        }
        coefficients.resize(m, false);
    }

    GF2mElement operator+(const GF2mElement& other) const {
        std::vector<bool> result_coeffs(m);
        for (int i = 0; i < m; ++i) {
            result_coeffs[i] = coefficients[i] ^ other.coefficients[i];
        }
        return GF2mElement(result_coeffs);
    }

    GF2mElement squareONB() const {
        std::vector<bool> squared_coeffs(m);
        squared_coeffs[m - 1] = coefficients[0];
        for (int i = 0; i < m - 1; ++i) {
            squared_coeffs[i] = coefficients[i + 1];
        }
        return GF2mElement(squared_coeffs);
    }

    bool trace() const {
        bool trace_value = false;
        for (bool coeff : coefficients) {
            trace_value ^= coeff;
        }
        return trace_value;
    }

    friend std::ostream& operator<<(std::ostream& os, const GF2mElement& element) {
        for (int i = m - 1; i >= 0; --i) {
            os << (element.coefficients[i] ? '1' : '0');
        }
        return os;
    }

    static int mod_pow_2(int exponent, int mod) {
        if (mod_pow_2_cache.find(exponent) != mod_pow_2_cache.end()) {
            return mod_pow_2_cache[exponent];
        }

        int result = 1;
        for (int i = 0; i < exponent; ++i) {
            result = (result * 2) % mod;
        }

        mod_pow_2_cache[exponent] = result;
        return result;
    }

    static int compute_matrix_element(int i, int j) {
        int mod = p;
        int pow_i = mod_pow_2(i, mod);
        int pow_j = mod_pow_2(j, mod);

        int results[] = {
                (pow_i + pow_j) % mod,
                (pow_i - pow_j + mod) % mod,
                (-pow_i + pow_j + mod) % mod,
                (-pow_i - pow_j + mod) % mod
        };

        for (int result : results) {
            if (result == 1 || result == -466) return 1;
        }

        return 0;
    }

    static std::vector<std::pair<int, int>> createMultiplicativeMatrix() {
        std::vector<std::pair<int, int>> one_positions;
        char prevVal = 0;
        for (int i = 0; i < m; ++i) {
            prevVal = 0;
            for (int j = 0; j < m; ++j) {
                if (compute_matrix_element(i, j) == 1) {
                    one_positions.emplace_back(i, j);
                    prevVal++;
                    if (prevVal == 2) break;
                }
            }
        }

        return one_positions;
    }

    std::vector<bool> transposeToVector() const {
        std::vector<bool> transposed_vector(m);
        for (int i = 0; i < m; ++i) {
            transposed_vector[i] = coefficients[m - 1 - i];
        }
        return transposed_vector;
    }

    static void printMatrix(const std::vector<std::vector<bool>>& matrix) {
        for (const auto& row : matrix) {
            for (bool val : row) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }

    std::vector<bool> multiplyWithMatrix(const std::vector<std::pair<int, int>>& one_positions) const {
        std::vector<bool> result(m, false);
        for (const auto& [i, j] : one_positions) {
            bool temp = result[i] ^ coefficients[m - 1 - j];
            result[i] = temp;
        }
        return result;
    }

    bool multiplyWithTransposed(const std::vector<bool>& other) const {
        bool result = false;
        for (int i = 0; i < m; ++i) {
            if (other[i]) {
                result ^= coefficients[i];
            }
        }
        return result;
    }

    GF2mElement cyclicLeftShift(int positions) const {
        int size = coefficients.size();
        std::vector<bool> shifted_coeffs(size);

        for (int i = 0; i < size; ++i) {
            int new_index = (i + positions) % size;
            shifted_coeffs[new_index] = coefficients[i];
        }

        return GF2mElement(shifted_coeffs);
    }

    void print() const {
        for (int i = m - 1; i >= 0; --i) {
            std::cout << (coefficients[i] ? '1' : '0');
        }
        std::cout << std::endl;
    }

    static std::string multiplyAndShift(GF2mElement a, GF2mElement b, int steps) {
        std::vector<std::pair<int, int>> matrixA = createMultiplicativeMatrix();
        std::string resultVector;

        for (int step = 0; step < steps; ++step) {
            std::vector<bool> productWithA = a.multiplyWithMatrix(matrixA);
            std::vector<bool> transposedB = b.transposeToVector();
            bool multiplicationResult = GF2mElement(productWithA).multiplyWithTransposed(transposedB);
            resultVector.push_back(multiplicationResult ? '1' : '0');

            a = a.cyclicLeftShift(1);
            b = b.cyclicLeftShift(1);
        }

        return resultVector;
    }

    GF2mElement operator*(const GF2mElement& other) const {
        std::string result = multiplyAndShift(*this, other, 233);
        std::vector<bool> result_vector;
        for (auto it = result.rbegin(); it != result.rend(); ++it) {
            result_vector.push_back(*it == '1');
        }
        return GF2mElement(result_vector);
    }

    GF2mElement power(const std::string& exponent) const {
        std::vector<bool> neutral_coeffs(m, true);
        GF2mElement result(neutral_coeffs);
        GF2mElement base = *this;

        if (!exponent.empty() && exponent[0] == '1') {
            result = result * base;
        }

        for (size_t i = 1; i < exponent.length(); ++i) {
            result = result.squareONB();
            if (exponent[i] == '1') {
                result = result * base;
            }

        }
        return result;
    }

    GF2mElement inverse() const {
        GF2mElement beta = *this;
        int k = 1;
        std::string m_binary = "11101000"; // m - 1= 232

        for (int i = 1; i <= 7; ++i) {
            GF2mElement original_beta = beta;
            for (int j = 0; j < k; ++j) {
                beta = beta.squareONB();
            }
            beta = beta * original_beta;
            k *= 2;

            if (m_binary[i] == '1') {
                GF2mElement squared_beta = beta.squareONB();
                beta = squared_beta * (*this);
                ++k;
            }
        }
        beta = beta.squareONB();
        return beta;
    }

};

std::unordered_map<int, int> GF2mElement::mod_pow_2_cache;

int main() {

    GF2mElement a("10111100000011111110110111100101101100100111011101101000001011110001001110001110110001011101100110100001001110101101011011100100000110011010111110010000001010100101111101010100000010011001001001110100110011101111100101011110010111010");
    GF2mElement b("10010100100111000111100100011001111101000111000010110011001110000101100000111110101110000100000001101110110001110001100101000111011010110111001101110111001000000101100101110011000001011010010101110111100111001010001000001111010001010");

//    std::cout << "a = " << a << std::endl;
//    std::cout << "b = " << b << std::endl << std::endl;
    std::cout << std::endl;

    std::string N = "00101001011111011010001010001101011000100101011011001110100011100111010111101101011000010111000100110011110011100100001001011101101110110101111111001010010001101011010100010010110001011001101100111111111011111100010010100011101000111";

    auto start_add = std::chrono::high_resolution_clock::now();
    GF2mElement c = a + b;
    auto stop_add = std::chrono::high_resolution_clock::now();
    auto duration_add = std::chrono::duration_cast<std::chrono::microseconds>(stop_add - start_add);
    std::cout << "Addition: " << c << std::endl;
    std::cout << "Time: " << duration_add.count() << " microseconds" << std::endl << std::endl;

    auto start_square = std::chrono::high_resolution_clock::now();
    GF2mElement a_squared = a.squareONB();
    auto stop_square = std::chrono::high_resolution_clock::now();
    auto duration_square = std::chrono::duration_cast<std::chrono::microseconds>(stop_square - start_square);
    std::cout << "a^2: " << a_squared << std::endl;
    std::cout << "Time: " << duration_square.count() << " microseconds" << std::endl << std::endl;

    auto start_trace = std::chrono::high_resolution_clock::now();
    bool trace_of_a = a.trace();
    auto stop_trace = std::chrono::high_resolution_clock::now();
    auto duration_trace = std::chrono::duration_cast<std::chrono::microseconds>(stop_trace - start_trace);
    std::cout << "Trace of a: " << trace_of_a << std::endl;
    std::cout << "Time: " << duration_trace.count() << " microseconds" << std::endl << std::endl;


    auto start_mul = std::chrono::high_resolution_clock::now();
    GF2mElement product = a * b;
    auto stop_mul = std::chrono::high_resolution_clock::now();
    auto duration_mul = std::chrono::duration_cast<std::chrono::microseconds>(stop_mul - start_mul);
    std::cout << "Multiplication: " << product << std::endl;
    std::cout << "Time: " << duration_mul.count() << " microseconds" << std::endl << std::endl;

    auto start_inv = std::chrono::high_resolution_clock::now();
    GF2mElement inverse_element = a.inverse();
    auto stop_inv = std::chrono::high_resolution_clock::now();
    auto duration_inv = std::chrono::duration_cast<std::chrono::microseconds>(stop_inv - start_inv);
    std::cout << "Inverse of a: " << inverse_element << std::endl;
    std::cout << "Time: " << duration_inv.count() << " microseconds" << std::endl << std::endl;

    auto start_pow = std::chrono::high_resolution_clock::now();
    GF2mElement a_pow = a.power(N);
    auto stop_pow = std::chrono::high_resolution_clock::now();
    auto duration_pow = std::chrono::duration_cast<std::chrono::microseconds>(stop_pow - start_pow);
    std::cout << "a^N : " << a_pow << std::endl;
    std::cout << "Time: " << duration_pow.count() << " microseconds" << std::endl;

    return 0;
}