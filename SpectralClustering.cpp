#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>

struct vec2 {
	float x, y;
	vec2(float x0 = 0, float y0 = 0) { x = x0; y = y0; }
	vec2 operator*(float a) const { return vec2(x * a, y * a); }
};

inline float dot(const vec2& v1, const vec2& v2) {
	return (v1.x * v2.x + v1.y * v2.y);
}

inline float length(const vec2& v) { return sqrtf(dot(v, v)); }
inline vec2 normalize(const vec2& v) { return v * (1 / length(v)); }

struct vec4 {
	float x, y, z, w;
	vec4(float x0 = 0, float y0 = 0, float z0 = 0, float w0 = 0) { x = x0; y = y0; z = z0; w = w0; }
};



class Matrix {
public:
	int x;
	int y;
	std::vector<std::vector<float>> data;
	Matrix(int sizex, int sizey) {
		x = sizex;
		y = sizey;
		data = std::vector<std::vector<float>>(x,std::vector<float>(y));
	}

	void identity() {
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < x; j++) {
				if (i == j) {
					data[i][j] = 1.0f;
				}
				else {
					data[i][j] = 0.0f;
				}
			}
		}
	}

	Matrix operator*(Matrix matrix) {
		Matrix m = Matrix(matrix.x, matrix.y);
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < x; j++) {
				m.data[i][j] = 0;
				for (int k = 0; k < x; k++) {
					m.data[i][j] += this->data[i][k] * matrix.data[k][j];
				}
				
			}
		}
		return m;
	}

	Matrix transposed(const Matrix& matrix) {
		Matrix m =  Matrix(matrix.x, matrix.y);
		for (int i = 0; i < matrix.x; i++) {
			for (int j = 0; j < matrix.y; j++) {
				m.data[i][j] = matrix.data[j][i];
			}
		}
		return m;
	}

};

Matrix eigenvectors2by2(float a, float b, float c, float d) {
	float e1 = (a + d) / 2.0f + sqrtf(b * b + ((a - d) / 2) * ((a - d) / 2));
	float e2 = (a + d) / 2.0f - sqrtf(b * b + ((a - d) / 2) * ((a - d) / 2));
	vec4 m2b2 = vec4(a - e1, b, c, d - e1);
	vec2 r = vec2(m2b2.y, -m2b2.x);
	r = normalize(r);
	Matrix m = Matrix(2, 2);
	m.data[0][0] = r.x;
	m.data[1][0] = r.y;
	m.data[0][1] = -r.y;
	m.data[1][1] = r.x;
	return m;
}

Matrix jacobialgoritm(Matrix m) {
	bool running = true;
	while (running) {
		int x = 0;
		int y = 1;
		float maxvalue = 0.0f;
		for (int i = 0; i < m.x; i++) {
			for (int j = 0; j < m.x; j++) {
				if (i != j) {
					if (maxvalue < abs(m.data[i][j])) {
						maxvalue = abs(m.data[i][j]);
						x = i;
						y = j;
					}
				}
			}
		}
		if (maxvalue < 0.00001f) {
			running = false;
		}
		else {

			Matrix u = eigenvectors2by2(m.data[x][x], m.data[x][y], m.data[y][x], m.data[y][y]);
			Matrix id = Matrix(m.x, m.y);
			id.identity();
			id.data[x][x] = u.data[0][0];
			id.data[x][y] = u.data[0][1];
			id.data[y][x] = u.data[1][0];
			id.data[y][y] = u.data[1][1];
			Matrix idt  = id.transposed(id);
			
			m = idt * m;
			m = m * id;
		
		}
	}
	return m;
}

std::vector<float> gaussELim(Matrix m, std::vector<float> b) {
	int n = b.size();
	std::vector<float> v(n);
	//elimination
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (m.data[j][i] != 0.0f) {
				float factor = m.data[i][i]/m.data[j][i];
				for (int k = i; k < n; k++) {
					m.data[j][k] = m.data[i][k] - m.data[j][k] * factor;
				}
				b[j] = b[i] - b[j] * factor;
			}
		}
	}

	//back-substitution
	if (m.data[n - 1][n - 1] != 0.0f) {
		v[n - 1] = b[n - 1] / m.data[n - 1][n - 1];
	}
	else {
		v[n - 1] = 1.0f;
	}
	
	for (int i = n - 2; i >= 0; i--) {
		float sum = 0.0f;
		for (int j = i + 1; j < n; j++) {
			sum = sum + m.data[i][j] * v[j];
		}
		if (m.data[i][i] == 0.0f) {
		//azt jeletni hogy lehet bármi?
			v[i] = 1.0f;
		}
		else {
			v[i] = (b[i] - sum) / m.data[i][i];
		}
	}
	return v;


}
int main(char* args) {
	Matrix m = Matrix(4, 4);
	m.data[0][0] = 2.0f;
	m.data[0][1] = 0.0f;
	m.data[0][2] = -1.0f;
	m.data[0][3] = -1.0f;
	m.data[1][0] = 0.0f;
	m.data[1][1] = 2.0f;
	m.data[1][2] = -1.0f;
	m.data[1][3] = -1.0f;
	m.data[2][0] = -1.0f;
	m.data[2][1] = -1.0f;
	m.data[2][2] = 2.0f;
	m.data[2][3] = 0.0f;
	m.data[3][0] = -1.0f;
	m.data[3][1] = -1.0f;
	m.data[3][2] = 0.0f;
	m.data[3][3] = 2.0f;

	m = jacobialgoritm(m);

	std::vector<float> v(4);
	for (int i = 0; i < 4; i++) {
		v[i] = 0.0f;
	}

	/*v = gaussELim(m, v);
	for (int i = 0; i < m.x; i++) {	
			printf("\t %f", v[i]);
		printf("\n");
	}*/
	for (int i = 0; i < m.x; i++) {
		for (int j = 0; j < m.x; j++) {
		printf("\t %f", m.data[i][j]);
		}
		printf("\n");
	}


	
}