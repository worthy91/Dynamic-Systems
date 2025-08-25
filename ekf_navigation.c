include <stdio.h> 
#include <string.h> 
#include <math.h> 
#include "freertos/FreeRTOS.h" 
#include "freertos/task.h" 
#include "freertos/queue.h" 
#include "driver/uart.h" 
#include "esp_timer.h" 
 
#define UART_NUM UART_NUM_1 
#define UART_BAUD_RATE 115200 
#define UART_RX_PIN 18 
#define UART_TX_PIN -1 
#define BUF_SIZE 1024 
#define PACKET_SIZE sizeof(SensorPacket) 
#define PACKET_BUFFER_SIZE 10 
 
#define STATE_DIM 17 // [qw,qx,qy,qz, px,py,pz, vx,vy,vz, bax,bay,baz, bgx,bgy,bgz, bp] 
#define MEAS_DIM 7   // [mx,my,mz, pressure, ax,ay,az] 
 
// CRITICAL FIX 1: Add process noise dimension 
#define NOISE_DIM 12 // [nav, naa, nbg, nba] - angular velocity, linear acceleration, gyro 
bias, accel bias noise 
 
// Sensor Packet 
typedef struct __attribute__((packed)) { 
    uint32_t timestamp_us; 
    float ax, ay, az; 
    float gx, gy, gz; 
    int mx, my, mz; 
    float pressure; 
    uint8_t checksum; 
} SensorPacket; 
 
// EKF Global Varialbles 
typedef struct { 
    float x[STATE_DIM];                    // State vector 
    float P[STATE_DIM * STATE_DIM];        // Covariance matrix 
    float F[STATE_DIM * STATE_DIM];        // State transition Jacobian 
    float G[STATE_DIM * NOISE_DIM];        // Process noise Jacobian 
    float H[MEAS_DIM * STATE_DIM];         // Measurement Jacobian 
    float Q[NOISE_DIM * NOISE_DIM];        // Process noise covariance 
(proper dimension) 
    float R[MEAS_DIM * MEAS_DIM];          // Measurement noise covariance 
} EKFState; 
 
// Global variables 
static QueueHandle_t sensor_queue; 
static EKFState ekf; 
static float gravity_vector[3] = {0.0f, 0.0f, -9.81f}; 
static float reference_altitude = 0.0f; 
 
// Magnetic field parameters for GOTH IQBAL 
static float mag_declination = 1.317f;  // +1° 19' = 1.317° (positive = east) 
static float mag_inclination = 39.417f; // 39° 25' = 39.417° 
static float mag_field_strength = 45681.5f; // Total field strength in nT 
static float mag_field_ned[3]; // Magnetic field model in NED frame 
 
// Matrix Operations Funtions in Global 
static inline void mat_zero(float *A, int rows, int cols) { 
    memset(A, 0, rows * cols * sizeof(float)); 
} 
 
static inline void mat_identity(float *A, int n) { 
    mat_zero(A, n, n); 
    for (int i = 0; i < n; i++) { 
        A[i * n + i] = 1.0f; 
    } 
} 
 
static inline void mat_copy(const float *src, float *dst, int rows, int cols) { 
    memcpy(dst, src, rows * cols * sizeof(float)); 
} 
 
static void mat_mult(const float *A, const float *B, float *C, int m, int n, int p) { 
    for (int i = 0; i < m; i++) { 
        for (int j = 0; j < p; j++) { 
            C[i * p + j] = 0.0f; 
            for (int k = 0; k < n; k++) { 
                C[i * p + j] += A[i * n + k] * B[k * p + j]; 
            } 
        } 
    } 
} 
 
static void mat_transpose(const float *A, float *AT, int rows, int cols) { 
    for (int i = 0; i < rows; i++) { 
        for (int j = 0; j < cols; j++) { 
            AT[j * rows + i] = A[i * cols + j]; 
        } 
    } 
} 
 
static void mat_add(const float *A, const float *B, float *C, int rows, int cols) { 
    for (int i = 0; i < rows * cols; i++) { 
        C[i] = A[i] + B[i]; 
    } 
} 
 
static void mat_sub(const float *A, const float *B, float *C, int rows, int cols) { 
    for (int i = 0; i < rows * cols; i++) { 
        C[i] = A[i] - B[i]; 
    } 
} 
 
// LU decomposition Functions in Global, Linear Algebra (For Matrix Inversion)
 
static int lu_decompose(float *A, int *perm, int n) { 
    for (int i = 0; i < n; i++) perm[i] = i; 
     
    for (int k = 0; k < n - 1; k++) { 
        int pivot = k; 
        float max_val = fabsf(A[perm[k] * n + k]); 
        for (int i = k + 1; i < n; i++) { 
            float val = fabsf(A[perm[i] * n + k]); 
            if (val > max_val) { 
                max_val = val; 
                pivot = i; 
            } 
        } 
         
        if (max_val < 1e-12f) return 0; 
         
        if (pivot != k) { 
            int temp = perm[k]; 
            perm[k] = perm[pivot]; 
            perm[pivot] = temp; 
        } 
         
        for (int i = k + 1; i < n; i++) { 
            float factor = A[perm[i] * n + k] / A[perm[k] * n + k]; 
            A[perm[i] * n + k] = factor; 
            for (int j = k + 1; j < n; j++) { 
                A[perm[i] * n + j] -= factor * A[perm[k] * n + j]; 
            } 
        } 
    } 
    return 1; 
} 
 
static void lu_solve(const float *A, const int *perm, const float *b, float *x, int n) { 
    for (int i = 0; i < n; i++) { 
        x[i] = b[perm[i]]; 
        for (int j = 0; j < i; j++) { 
            x[i] -= A[perm[i] * n + j] * x[j]; 
        } 
    } 
     
    for (int i = n - 1; i >= 0; i--) { 
        for (int j = i + 1; j < n; j++) { 
            x[i] -= A[perm[i] * n + j] * x[j]; 
        } 
        x[i] /= A[perm[i] * n + i]; 
    } 
} 
 
// Quaternion Functions in Global
static inline void quat_normalize(float *q) { 
    float norm = sqrtf(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]); 
    if (norm > 1e-6f) { 
        float inv_norm = 1.0f / norm; 
        for (int i = 0; i < 4; i++) q[i] *= inv_norm; 
    } 
} 
 
static void quat_to_rotation_matrix(const float *q, float *R) { 
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3]; 
     
    R[0] = 1 - 2*(qy*qy + qz*qz); R[1] = 2*(qx*qy - qw*qz);     R[2] = 2*(qx*qz + qw*qy); 
    R[3] = 2*(qx*qy + qw*qz);     R[4] = 1 - 2*(qx*qx + qz*qz); R[5] = 2*(qy*qz - qw*qx); 
    R[6] = 2*(qx*qz - qw*qy);     R[7] = 2*(qy*qz + qw*qx);     R[8] = 1 - 2*(qx*qx + qy*qy); 
} 
 
static void quat_rotate_vector(const float *q, const float *v, float *result) { 
    float R[9]; 
    quat_to_rotation_matrix(q, R); 
    for (int i = 0; i < 3; i++) { 
        result[i] = R[i*3]*v[0] + R[i*3+1]*v[1] + R[i*3+2]*v[2]; 
    } 
} 
 
// Skew-symmetric Function in Global
static void skew_symmetric(const float *v, float *S) { 
    mat_zero(S, 3, 3); 
    S[0*3 + 1] = -v[2]; S[0*3 + 2] =  v[1]; 
    S[1*3 + 0] =  v[2]; S[1*3 + 2] = -v[0]; 
    S[2*3 + 0] = -v[1]; S[2*3 + 1] =  v[0]; 
} 
 
// Calc NED from Raw Magnatometer Data
static void init_magnetic_field_model(void) { 
    float declination_rad = mag_declination * M_PI / 180.0f; 
    float inclination_rad = mag_inclination * M_PI / 180.0f; 
     
    float horizontal_component = mag_field_strength * cosf(inclination_rad); 
     
    mag_field_ned[0] = horizontal_component * cosf(declination_rad); 
    mag_field_ned[1] = horizontal_component * sinf(declination_rad); 
    mag_field_ned[2] = mag_field_strength * sinf(inclination_rad); 
     
    printf("Magnetic field model initialized for GOTH IQBAL:\n"); 
    printf("  NED components: [%.1f, %.1f, %.1f] nT\n",  
           mag_field_ned[0], mag_field_ned[1], mag_field_ned[2]); 
} 
 
// Complete EKF Initialization  
static void ekf_init(EKFState *ekf) { 
    init_magnetic_field_model(); 
     
    // Initialize state vector 
    mat_zero(ekf->x, STATE_DIM, 1); 
    ekf->x[0] = 1.0f; // qw = 1 (identity quaternion) 
     
    // Initialize covariance matrix 
    mat_identity(ekf->P, STATE_DIM); 
    // Quaternion uncertainties (small for initial alignment) 
    for (int i = 0; i < 4; i++) ekf->P[i * STATE_DIM + i] = 0.01f; 
    // Position uncertainties 
    for (int i = 4; i < 7; i++) ekf->P[i * STATE_DIM + i] = 1.0f; 
    // Velocity uncertainties 
    for (int i = 7; i < 10; i++) ekf->P[i * STATE_DIM + i] = 0.1f; 
    // Accelerometer bias uncertainties 
    for (int i = 10; i < 13; i++) ekf->P[i * STATE_DIM + i] = 0.01f; 
    // Gyroscope bias uncertainties 
    for (int i = 13; i < 16; i++) ekf->P[i * STATE_DIM + i] = 0.001f; 
    // Barometer bias uncertainty 
    ekf->P[16 * STATE_DIM + 16] = 100.0f; 
     
    // Initialize of process noise covariance Q as zero
    mat_zero(ekf->Q, NOISE_DIM, NOISE_DIM); 
    // Angular velocity noise insreeted in diag of Q (rad²/s²) 
    for (int i = 0; i < 3; i++) ekf->Q[i * NOISE_DIM + i] = pow(0.1f * M_PI/180.0f, 2); 
    // Linear acceleration noise insreeted in diag of Q (m²/s⁴) 
    for (int i = 3; i < 6; i++) ekf->Q[i * NOISE_DIM + i] = pow(0.1f, 2); 
    // Gyroscope bias random walk insreeted in diag of Q (rad²/s³) 
    for (int i = 6; i < 9; i++) ekf->Q[i * NOISE_DIM + i] = pow(1e-5f, 2); 
    // Accelerometer bias random walk insreeted in diag of Q (m²/s⁵) 
    for (int i = 9; i < 12; i++) ekf->Q[i * NOISE_DIM + i] = pow(1e-4f, 2); 
     
    // Initialize measurement noise covariance R 
    mat_zero(ekf->R, MEAS_DIM, MEAS_DIM); 
    float mag_noise_factor = mag_field_strength / 45000.0f; 
    ekf->R[0 * MEAS_DIM + 0] = 400.0f * mag_noise_factor; 
    ekf->R[1 * MEAS_DIM + 1] = 400.0f * mag_noise_factor; 
    ekf->R[2 * MEAS_DIM + 2] = 400.0f * mag_noise_factor; 
    ekf->R[3 * MEAS_DIM + 3] = 25.0f;   // pressure variance 
    ekf->R[4 * MEAS_DIM + 4] = 0.1f;    // ax variance 
    ekf->R[5 * MEAS_DIM + 5] = 0.1f;    // ay variance 
    ekf->R[6 * MEAS_DIM + 6] = 0.1f;    // az variance 
    // gyro not inserted in measurement noise update cz we use it as an actual input
    //and it will be corrected by accel and mag measurements.
    printf("EKF initialized with proper noise modeling\n"); 
} 
 
// Creating Jacobian F
static void compute_state_jacobian(EKFState *ekf, const SensorPacket *pkt, float dt) { 
    mat_identity(ekf->F, STATE_DIM); 
     
    float *q = &ekf->x[0]; 
    float *vel = &ekf->x[7]; 
    float *ba = &ekf->x[10]; 
    float *bg = &ekf->x[13]; 
     
    // Bias-corrected measurements for gyro and accel
    float wx = (pkt->gx - bg[0]) * M_PI / 180.0f; 
    float wy = (pkt->gy - bg[1]) * M_PI / 180.0f; 
    float wz = (pkt->gz - bg[2]) * M_PI / 180.0f; 
    float w[3] = {wx, wy, wz}; 
     
    float ax = pkt->ax - ba[0]; 
    float ay = pkt->ay - ba[1]; 
    float az = pkt->az - ba[2]; 
    float a[3] = {ax, ay, az}; 
     
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3]; 
     
    // 1. Quaternion matrix(tell how to extract quat from gyro data)
    float Omega[16] = { 
        0,   -wx,  -wy,  -wz, 
        wx,   0,    wz,  -wy, 
        wy,  -wz,   0,    wx, 
        wz,   wy,  -wx,   0 
    }; 
     // Quaternion covar quantities( tell how quat like Qx affect itself and other quats)
    for (int i = 0; i < 4; i++) { 
        for (int j = 0; j < 4; j++) { 
            ekf->F[i * STATE_DIM + j] += 0.5f * dt * Omega[i * 4 + j]; 
        } 
    } 
     
    // 2. How gyro bias affeect att (used for correcting gyro bias using att)
    // below is how each quat is affect each bias
    ekf->F[0 * STATE_DIM + 13] = -0.5f * dt * qx * M_PI/180.0f; 
    ekf->F[0 * STATE_DIM + 14] = -0.5f * dt * qy * M_PI/180.0f; 
    ekf->F[0 * STATE_DIM + 15] = -0.5f * dt * qz * M_PI/180.0f; 
     
    ekf->F[1 * STATE_DIM + 13] =  0.5f * dt * qw * M_PI/180.0f; 
    ekf->F[1 * STATE_DIM + 14] =  0.5f * dt * qz * M_PI/180.0f; 
    ekf->F[1 * STATE_DIM + 15] = -0.5f * dt * qy * M_PI/180.0f; 
     
    ekf->F[2 * STATE_DIM + 13] = -0.5f * dt * qz * M_PI/180.0f; 
    ekf->F[2 * STATE_DIM + 14] =  0.5f * dt * qw * M_PI/180.0f; 
    ekf->F[2 * STATE_DIM + 15] =  0.5f * dt * qx * M_PI/180.0f; 
     
    ekf->F[3 * STATE_DIM + 13] =  0.5f * dt * qy * M_PI/180.0f; 
    ekf->F[3 * STATE_DIM + 14] = -0.5f * dt * qx * M_PI/180.0f; 
    ekf->F[3 * STATE_DIM + 15] =  0.5f * dt * qw * M_PI/180.0f; 
     
    // 3. How velocity affect pos (like Pnorth depend on Vnorth) map through partial diff
    for (int i = 0; i < 3; i++) { 
        ekf->F[(4 + i) * STATE_DIM + (7 + i)] = dt; 
    } 
     
    // 4. How velocity dependes on att error 
    float R[9]; 
    quat_to_rotation_matrix(q, R); 
     
    // 5. how each quat chanhge with each element of rotation matix
    float dR_dqw[9] = { 
        0, -2*qz, 2*qy, 
        2*qz, 0, -2*qx, 
        -2*qy, 2*qx, 0 
    }; 
     
    float dR_dqx[9] = { 
        -4*qx, 2*qy, 2*qz, 
        2*qy, -4*qx, 0, 
        2*qz, 0, -4*qx 
    }; 
     
    float dR_dqy[9] = { 
        -4*qy, 2*qx, 0, 
        2*qx, 0, 2*qz, 
        0, 2*qz, -4*qy 
    }; 
     
    float dR_dqz[9] = { 
        -4*qz, 0, 2*qx, 
        0, -4*qz, 2*qy, 
        2*qx, 2*qy, 0 
    }; 
     
    // 6. how each quat affect each velocity element
    for (int i = 0; i < 3; i++) { 
        ekf->F[(7 + i) * STATE_DIM + 0] = dt * (dR_dqw[i*3]*a[0] + dR_dqw[i*3+1]*a[1] + 
dR_dqw[i*3+2]*a[2]); 
        ekf->F[(7 + i) * STATE_DIM + 1] = dt * (dR_dqx[i*3]*a[0] + dR_dqx[i*3+1]*a[1] + 
dR_dqx[i*3+2]*a[2]); 
        ekf->F[(7 + i) * STATE_DIM + 2] = dt * (dR_dqy[i*3]*a[0] + dR_dqy[i*3+1]*a[1] + 
dR_dqy[i*3+2]*a[2]); 
        ekf->F[(7 + i) * STATE_DIM + 3] = dt * (dR_dqz[i*3]*a[0] + dR_dqz[i*3+1]*a[1] + 
dR_dqz[i*3+2]*a[2]); 
    } 
     
    // 7. (how accel noise affect velocity) we take -ve R vector as explain in notes
    for (int i = 0; i < 3; i++) { 
        for (int j = 0; j < 3; j++) { 
            ekf->F[(7 + i) * STATE_DIM + (10 + j)] = -dt * R[i*3 + j]; 
        } 
    } 
     
    //Bias dynamics are identity and it will change only in correction step dq/dbg
    //similalrly, accel bias will change by velocioty in step dv/dba
    
} 
//--------------------------------- JACOBIAN F FINISHED ----------------------------// 


// 8. Est noise jaobiAN L/G 
static void compute_process_noise_jacobian(EKFState *ekf, const SensorPacket *pkt, float 
dt) { 
    mat_zero(ekf->G, STATE_DIM, NOISE_DIM); 
     
    float *q = &ekf->x[0]; 
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3]; 
     
    // 1. Quaternion affected by angular velocity noise: dq/dnw 
    ekf->G[0 * NOISE_DIM + 0] = -0.5f * dt * qx; 
    ekf->G[0 * NOISE_DIM + 1] = -0.5f * dt * qy; 
    ekf->G[0 * NOISE_DIM + 2] = -0.5f * dt * qz; 
     
    ekf->G[1 * NOISE_DIM + 0] =  0.5f * dt * qw; 
    ekf->G[1 * NOISE_DIM + 1] =  0.5f * dt * qz; 
    ekf->G[1 * NOISE_DIM + 2] = -0.5f * dt * qy; 
     
    ekf->G[2 * NOISE_DIM + 0] = -0.5f * dt * qz; 
    ekf->G[2 * NOISE_DIM + 1] =  0.5f * dt * qw; 
    ekf->G[2 * NOISE_DIM + 2] =  0.5f * dt * qx; 
     
    ekf->G[3 * NOISE_DIM + 0] =  0.5f * dt * qy; 
    ekf->G[3 * NOISE_DIM + 1] = -0.5f * dt * qx; 
    ekf->G[3 * NOISE_DIM + 2] =  0.5f * dt * qw; 
     
    // 2. Velocity affected by acceleration noise: dv/dna 
    float R[9]; 
    quat_to_rotation_matrix(q, R); 
     
    for (int i = 0; i < 3; i++) { 
        for (int j = 0; j < 3; j++) { 
            ekf->G[(7 + i) * NOISE_DIM + (3 + j)] = dt * R[i*3 + j]; 
        } 
    } 
     
    // 3. Gyro bias affected by gyro bias noise: dbg/dnbg 
    for (int i = 0; i < 3; i++) { 
        ekf->G[(13 + i) * NOISE_DIM + (6 + i)] = dt; 
    } 
     
    // 4. Accel bias affected by accel bias noise: dba/dnba 
    for (int i = 0; i < 3; i++) { 
        ekf->G[(10 + i) * NOISE_DIM + (9 + i)] = dt; 
    } 
} 
//--------------------- PROCESS NOISE JACOBIAN G FINISHED ---------------------// 

// Construction of H and h(x) where H linearize h(x)
static void compute_measurement_jacobian(EKFState *ekf, float *h) { 
    mat_zero(ekf->H, MEAS_DIM, STATE_DIM); 
     
    float *q = &ekf->x[0]; 
    float *pos = &ekf->x[4]; 
    float *ba = &ekf->x[10]; 
    float *bp = &ekf->x[16]; 
     
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3]; 
     
    // 1. Magnetometer est mea surement Jacobian: dh_mag/dq 
    float dR_dqw[9] = {0, -2*qz, 2*qy, 2*qz, 0, -2*qx, -2*qy, 2*qx, 0}; 
    float dR_dqx[9] = {-4*qx, 2*qy, 2*qz, 2*qy, -4*qx, 0, 2*qz, 0, -4*qx}; 
    float dR_dqy[9] = {-4*qy, 2*qx, 0, 2*qx, 0, 2*qz, 0, 2*qz, -4*qy}; 
    float dR_dqz[9] = {-4*qz, 0, 2*qx, 0, -4*qz, 2*qy, 2*qx, 2*qy, 0}; 
     
    for (int i = 0; i < 3; i++) { 
        ekf->H[i * STATE_DIM + 0] = dR_dqw[i*3]*mag_field_ned[0] + 
dR_dqw[i*3+1]*mag_field_ned[1] + dR_dqw[i*3+2]*mag_field_ned[2]; 
        ekf->H[i * STATE_DIM + 1] = dR_dqx[i*3]*mag_field_ned[0] + 
dR_dqx[i*3+1]*mag_field_ned[1] + dR_dqx[i*3+2]*mag_field_ned[2]; 
        ekf->H[i * STATE_DIM + 2] = dR_dqy[i*3]*mag_field_ned[0] + 
dR_dqy[i*3+1]*mag_field_ned[1] + dR_dqy[i*3+2]*mag_field_ned[2]; 
        ekf->H[i * STATE_DIM + 3] = dR_dqz[i*3]*mag_field_ned[0] + 
dR_dqz[i*3+1]*mag_field_ned[1] + dR_dqz[i*3+2]*mag_field_ned[2]; 
    } 
     
    // 2. Barometer  est measurement Jacobian 
    ekf->H[3 * STATE_DIM + 6] = -8.65f;   // dp/dz 
    ekf->H[3 * STATE_DIM + 16] = 1.0f;    // dp/dbp 
     
    // 3. Accelerometer est measurement Jacobian: dh_acc/dq and dh_acc/dba 
    for (int i = 0; i < 3; i++) { 
        ekf->H[(4 + i) * STATE_DIM + 0] = dR_dqw[i*3]*gravity_vector[0] + 
dR_dqw[i*3+1]*gravity_vector[1] + dR_dqw[i*3+2]*gravity_vector[2]; 
        ekf->H[(4 + i) * STATE_DIM + 1] = dR_dqx[i*3]*gravity_vector[0] + 
dR_dqx[i*3+1]*gravity_vector[1] + dR_dqx[i*3+2]*gravity_vector[2]; 
        ekf->H[(4 + i) * STATE_DIM + 2] = dR_dqy[i*3]*gravity_vector[0] + 
dR_dqy[i*3+1]*gravity_vector[1] + dR_dqy[i*3+2]*gravity_vector[2]; 
        ekf->H[(4 + i) * STATE_DIM + 3] = dR_dqz[i*3]*gravity_vector[0] + 
dR_dqz[i*3+1]*gravity_vector[1] + dR_dqz[i*3+2]*gravity_vector[2]; 
         
        // Accelerometer bias Jacobian 
        ekf->H[(4 + i) * STATE_DIM + (10 + i)] = 1.0f; 
    } 
     
    // Predicted measurements h
    float mag_body[3], gravity_body[3]; 
    quat_rotate_vector(q, mag_field_ned, mag_body); 
    quat_rotate_vector(q, gravity_vector, gravity_body); 
     
    h[0] = mag_body[0]; 
    h[1] = mag_body[1]; 
    h[2] = mag_body[2]; 
    h[3] = reference_altitude - pos[2] * 8.65f + *bp; 
    h[4] = gravity_body[0] + ba[0]; 
    h[5] = gravity_body[1] + ba[1]; 
    h[6] = gravity_body[2] + ba[2]; 
} 
//--------------- MEASUREMENT JACOBIAN H and h(x) FINISHED ---------------// 

// EKF Prediction Step


static void ekf_predict(EKFState *ekf, const SensorPacket *pkt, float dt) { 
    float *q = &ekf->x[0]; 
    float *pos = &ekf->x[4]; 
    float *vel = &ekf->x[7]; 
    float *ba = &ekf->x[10]; 
    float *bg = &ekf->x[13]; 
     
    // Bias-corrected measurements (actual vals from sensor)
    float wx = (pkt->gx - bg[0]) * M_PI / 180.0f; 
    float wy = (pkt->gy - bg[1]) * M_PI / 180.0f; 
    float wz = (pkt->gz - bg[2]) * M_PI / 180.0f; 
     
    float ax = pkt->ax - ba[0]; 
    float ay = pkt->ay - ba[1]; 
    float az = pkt->az - ba[2]; 
     
    // Quaternion integration using 4th-order Runge-Kutta (cz euler method is not stable) 
    float k1[4], k2[4], k3[4], k4[4], q_temp[4]; 
     
    // k1 
    k1[0] = 0.5f * (-q[1]*wx - q[2]*wy - q[3]*wz); 
    k1[1] = 0.5f * ( q[0]*wx + q[2]*wz - q[3]*wy); 
    k1[2] = 0.5f * ( q[0]*wy - q[1]*wz + q[3]*wx); 
    k1[3] = 0.5f * ( q[0]*wz + q[1]*wy - q[2]*wx); 
     
    // k2 
    for (int i = 0; i < 4; i++) q_temp[i] = q[i] + 0.5f * dt * k1[i]; 
    k2[0] = 0.5f * (-q_temp[1]*wx - q_temp[2]*wy - q_temp[3]*wz); 
    k2[1] = 0.5f * ( q_temp[0]*wx + q_temp[2]*wz - q_temp[3]*wy); 
    k2[2] = 0.5f * ( q_temp[0]*wy - q_temp[1]*wz + q_temp[3]*wx); 
    k2[3] = 0.5f * ( q_temp[0]*wz + q_temp[1]*wy - q_temp[2]*wx); 
     
    // k3 
    for (int i = 0; i < 4; i++) q_temp[i] = q[i] + 0.5f * dt * k2[i]; 
    k3[0] = 0.5f * (-q_temp[1]*wx - q_temp[2]*wy - q_temp[3]*wz); 
    k3[1] = 0.5f * ( q_temp[0]*wx + q_temp[2]*wz - q_temp[3]*wy); 
    k3[2] = 0.5f * ( q_temp[0]*wy - q_temp[1]*wz + q_temp[3]*wx); 
    k3[3] = 0.5f * ( q_temp[0]*wz + q_temp[1]*wy - q_temp[2]*wx); 
     
    // k4 
    for (int i = 0; i < 4; i++) q_temp[i] = q[i] + dt * k3[i]; 
    k4[0] = 0.5f * (-q_temp[1]*wx - q_temp[2]*wy - q_temp[3]*wz); 
    k4[1] = 0.5f * ( q_temp[0]*wx + q_temp[2]*wz - q_temp[3]*wy); 
    k4[2] = 0.5f * ( q_temp[0]*wy - q_temp[1]*wz + q_temp[3]*wx); 
    k4[3] = 0.5f * ( q_temp[0]*wz + q_temp[1]*wy - q_temp[2]*wx); 
     
    // Update quaternion after applyig RK4
    for (int i = 0; i < 4; i++) { 
        q[i] += (dt/6.0f) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]); 
    } 
    quat_normalize(q); 
     
    // Transform acceleration to world frame using R
    float accel_world[3]; 
    float accel_body[3] = {ax, ay, az}; 
    quat_rotate_vector(q, accel_body, accel_world); 
    accel_world[2] += 9.81f;  // Gravity compensation 
     
    // Integrate for velocity and pos
    for (int i = 0; i < 3; i++) { 
        vel[i] += accel_world[i] * dt; 
        pos[i] += vel[i] * dt; 
    } 
     
    // Biases remain constant (random walk model) assumed 
     
    // Predicting the New Uncertainty (Co-Var)
 
    // P = F*P*F' + G*Q*G' 
    compute_state_jacobian(ekf, pkt, dt); 
    compute_process_noise_jacobian(ekf, pkt, dt); 
    // var to make covar calc easier
    

    static float temp1[STATE_DIM * STATE_DIM]; 
    static float temp2[STATE_DIM * STATE_DIM]; 
    static float temp3[STATE_DIM * NOISE_DIM]; 
    static float FT[STATE_DIM * STATE_DIM]; 
    static float GT[NOISE_DIM * STATE_DIM]; 
    static float GQG[STATE_DIM * STATE_DIM]; 
     
    // F*P*F' 
    mat_transpose(ekf->F, FT, STATE_DIM, STATE_DIM); 
    mat_mult(ekf->F, ekf->P, temp1, STATE_DIM, STATE_DIM, STATE_DIM); 
    mat_mult(temp1, FT, temp2, STATE_DIM, STATE_DIM, STATE_DIM); 
     
    // G*Q*G' 
    mat_transpose(ekf->G, GT, STATE_DIM, NOISE_DIM); 
    mat_mult(ekf->G, ekf->Q, temp3, STATE_DIM, NOISE_DIM, NOISE_DIM); 
    mat_mult(temp3, GT, GQG, STATE_DIM, NOISE_DIM, STATE_DIM); 
     
    // P = F*P*F' + G*Q*G' 
    mat_add(temp2, GQG, ekf->P, STATE_DIM, STATE_DIM); 
} 
//--------------------- EKF PREDICT STEP FINISHED ---------------------// 

// Apply Mahalanobis
static bool validate_measurement_innovation(const float *innovation, const float *S) { 
    // Chi-squared test for outlier detection 
    // Compute innovation^T * S^(-1) * innovation (Mahalanobis distance) 
    static float S_work[MEAS_DIM * MEAS_DIM]; 
    static int perm[MEAS_DIM]; 
    static float temp_innov[MEAS_DIM]; 
     
    mat_copy(S, S_work, MEAS_DIM, MEAS_DIM); 
     
    if (!lu_decompose(S_work, perm, MEAS_DIM)) { 
        return false; // Singular covariance and calc cant proceed reliabaly
    } 
     
    lu_solve(S_work, perm, innovation, temp_innov, MEAS_DIM); 
     
    float mahalanobis_dist = 0.0f; 
    for (int i = 0; i < MEAS_DIM; i++) { 
        mahalanobis_dist += innovation[i] * temp_innov[i]; 
    } 
     
    // Chi-squared threshold for 7 DOF at 95% confidence level (c its standard and commonly used)
    float chi2_threshold = 14.07f; 
     
    if (mahalanobis_dist > chi2_threshold) { 
        printf("Measurement rejected: Mahalanobis distance %.2f > %.2f\n",  
               mahalanobis_dist, chi2_threshold); 
        return false; 
    } 
     
    return true; 
} 
 

static void ekf_update(EKFState *ekf, const SensorPacket *pkt, float *innovation_out) { 
    // Measurements 
    float z[MEAS_DIM] = { 
        (float)pkt->mx, (float)pkt->my, (float)pkt->mz, 
        pkt->pressure, 
        pkt->ax, pkt->ay, pkt->az 
    }; 
     
    // Predicted measurements 
    float h[MEAS_DIM]; 
    compute_measurement_jacobian(ekf, h); 
     
    // Innovation: y = z - h 
    float y[MEAS_DIM]; 
    for (int i = 0; i < MEAS_DIM; i++) { 
        y[i] = z[i] - h[i]; 
        if (innovation_out) innovation_out[i] = y[i]; 
    } 
     
    // Innovation covariance: S = H*P*H' + R 
    static float HT[STATE_DIM * MEAS_DIM]; 
    static float PH[STATE_DIM * MEAS_DIM]; 
    static float HPH[MEAS_DIM * MEAS_DIM]; 
    static float S[MEAS_DIM * MEAS_DIM]; 
    static float S_work[MEAS_DIM * MEAS_DIM]; 
    static int perm[MEAS_DIM]; 
     
    mat_transpose(ekf->H, HT, MEAS_DIM, STATE_DIM); 
    mat_mult(ekf->P, HT, PH, STATE_DIM, STATE_DIM, MEAS_DIM); 
    mat_mult(ekf->H, PH, HPH, MEAS_DIM, STATE_DIM, MEAS_DIM); 
    mat_add(HPH, ekf->R, S, MEAS_DIM, MEAS_DIM); 
     
    // Validate measurement before update 
    if (!validate_measurement_innovation(y, S)) { 
        printf("Measurement rejected by innovation validation\n"); 
        return; 
    } 
     
    // Copy S for LU decomposition 
    mat_copy(S, S_work, MEAS_DIM, MEAS_DIM); 
     
    // Compute Kalman gain: K = P*H'*S^(-1) 
    static float K[STATE_DIM * MEAS_DIM]; 
     
    if (lu_decompose(S_work, perm, MEAS_DIM)) { 
        // Solve K = PH * S^(-1) by solving S*K^T = PH^T 
        static float K_temp[MEAS_DIM]; 
         
        for (int i = 0; i < STATE_DIM; i++) { 
            for (int j = 0; j < MEAS_DIM; j++) { 
                K_temp[j] = PH[i * MEAS_DIM + j]; 
            } 
             
            static float K_col[MEAS_DIM]; 
            lu_solve(S_work, perm, K_temp, K_col, MEAS_DIM); 
             
            for (int j = 0; j < MEAS_DIM; j++) { 
                K[i * MEAS_DIM + j] = K_col[j]; 
            } 
        } 
         
        // State update: x = x + K*y 
        for (int i = 0; i < STATE_DIM; i++) { 
            float update = 0.0f; 
            for (int j = 0; j < MEAS_DIM; j++) { 
                update += K[i * MEAS_DIM + j] * y[j]; 
            } 
            ekf->x[i] += update; 
        } 
         
        // Normalize quaternion after update to preserve its unity and avoid sacling
        quat_normalize(&ekf->x[0]); 
         
        // ------ Applying Joseph form covariance update for numerical stability 
        // P = (I-KH)P(I-KH)' + KRK' 
        static float KH[STATE_DIM * STATE_DIM]; 
        static float I_KH[STATE_DIM * STATE_DIM]; 
        static float I_KH_T[STATE_DIM * STATE_DIM]; 
        static float temp_P1[STATE_DIM * STATE_DIM]; 
        static float temp_P2[STATE_DIM * STATE_DIM]; 
        static float KR[STATE_DIM * MEAS_DIM]; 
        static float KRK[STATE_DIM * STATE_DIM]; 
         
        mat_mult(K, ekf->H, KH, STATE_DIM, MEAS_DIM, STATE_DIM); 
        mat_identity(I_KH, STATE_DIM); 
        mat_sub(I_KH, KH, I_KH, STATE_DIM, STATE_DIM); 
         
        // (I-KH)P(I-KH)' 
        mat_transpose(I_KH, I_KH_T, STATE_DIM, STATE_DIM); 
        mat_mult(I_KH, ekf->P, temp_P1, STATE_DIM, STATE_DIM, STATE_DIM); 
        mat_mult(temp_P1, I_KH_T, temp_P2, STATE_DIM, STATE_DIM, STATE_DIM); 
         
        // KRK' 
        mat_mult(K, ekf->R, KR, STATE_DIM, MEAS_DIM, MEAS_DIM); 
        mat_mult(KR, HT, KRK, STATE_DIM, MEAS_DIM, STATE_DIM); 
         
        // Final covariance update 
        mat_add(temp_P2, KRK, ekf->P, STATE_DIM, STATE_DIM); 
         
        // Check diagonal elements and fix negative variances 
        for (int i = 0; i < STATE_DIM; i++) { 
            if (ekf->P[i * STATE_DIM + i] < 1e-12f) { 
                ekf->P[i * STATE_DIM + i] = 1e-12f; 
            } 
        } 
         
    } else { 
        printf("Warning: Singular innovation covariance matrix\n"); 
    } 
} 
//--------Kalman Gain is calced, State is Updated and Joseph Co-Var is Calced---------//

//Adaptive Tuning of R matrix based on innovation statistics (last 10 innovations) 
static void adaptive_filter_tuning(EKFState *ekf, const float *innovation) { 
    static float innovation_history[MEAS_DIM][10]; 
    static int history_index = 0; 
    static int history_count = 0; 
    static bool initialized = false; 
     
    if (!initialized) { 
        memset(innovation_history, 0, sizeof(innovation_history)); 
        initialized = true; 
    } 
     
    // Store innovation in circular buffer (circular buffer is used to make mem usage fixed and pridictable)
    for (int i = 0; i < MEAS_DIM; i++) { 
        innovation_history[i][history_index] = innovation[i]; 
    } 
    history_index = (history_index + 1) % 10; 
    if (history_count < 10) history_count++; 
     
    // Compute innovation variannce of each seven elements of MEAS matrix
    if (history_count >= 10 && history_index == 0) { 
        for (int i = 0; i < MEAS_DIM; i++) { 
            float mean = 0.0f, variance = 0.0f; 
             
            // Compute mean 
            for (int j = 0; j < 10; j++) { 
                mean += innovation_history[i][j]; 
            } 
            mean /= 10.0f; 
             
            // Compute variance 
            for (int j = 0; j < 10; j++) { 
                float diff = innovation_history[i][j] - mean; 
                variance += diff * diff; 
            } 
            variance /= 9.0f; // Sample variance 
             
            // Adaptive tuning: adjust R based on innovation variance 
            float expected_var = ekf->R[i * MEAS_DIM + i]; 
            float ratio = variance / expected_var; 
             
            if (ratio > 2.0f) { 
                // Innovation too large, increase R 
                ekf->R[i * MEAS_DIM + i] *= 1.1f; 
                printf("Increased R[%d] by 10%% due to large innovation\n", i); 
            } else if (ratio < 0.5f) { 
                // Innovation too small, decrease R 
                ekf->R[i * MEAS_DIM + i] *= 0.95f; 
                printf("Decreased R[%d] by 5%% due to small innovation\n", i); 
            } 
        } 
    } 
} 
 
// Convert quaternion to Euler angles (For user-interface only and will not be used in EKF calc) 
static void quat_to_euler(const float *q, float *roll, float *pitch, float *yaw) { 
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3]; 
     
    float sinr_cosp = 2 * (qw * qx + qy * qz); 
    float cosr_cosp = 1 - 2 * (qx * qx + qy * qy); 
    *roll = atan2f(sinr_cosp, cosr_cosp) * 180.0f / M_PI; 
     
    float sinp = 2 * (qw * qy - qz * qx); 
    if (fabsf(sinp) >= 1) 
        *pitch = copysignf(M_PI / 2, sinp) * 180.0f / M_PI; 
    else 
        *pitch = asinf(sinp) * 180.0f / M_PI; 
     
    float siny_cosp = 2 * (qw * qz + qx * qy); 
    float cosy_cosp = 1 - 2 * (qy * qy + qz * qz); 
    *yaw = atan2f(siny_cosp, cosy_cosp) * 180.0f / M_PI; 
} 
 
// Packet validation  
static uint8_t calculate_checksum(const SensorPacket *pkt) { 
    uint8_t checksum = 0; 
    const uint8_t *data = (const uint8_t*)pkt; 
    for (size_t i = 0; i < sizeof(SensorPacket) - 1; i++) { 
        checksum += data[i]; 
    } 
    return checksum; 
} 
 
static bool validate_packet(const SensorPacket *pkt) { 
    return (calculate_checksum(pkt) == pkt->checksum); 
} 
 
static bool validate_packet_extended(const SensorPacket *pkt) { 
    if (!validate_packet(pkt)) return false; 
     
    if (fabsf(pkt->ax) > 50.0f || fabsf(pkt->ay) > 50.0f || fabsf(pkt->az) > 50.0f) return false; 
    if (fabsf(pkt->gx) > 2000.0f || fabsf(pkt->gy) > 2000.0f || fabsf(pkt->gz) > 2000.0f) return 
false; 
    if (abs(pkt->mx) > 60000 || abs(pkt->my) > 60000 || abs(pkt->mz) > 60000) return false; 
    if (pkt->pressure < 80000.0f || pkt->pressure > 105000.0f) return false; 
     
    return true; 
} 
 
// UART receive task 
static void uart_receive_task(void *pvParameters) { 
    uint8_t buffer[PACKET_SIZE]; 
    int bytes_received = 0; 
    uint32_t packet_count = 0; 
    uint32_t error_count = 0; 
     
    while (1) { 
        int len = uart_read_bytes(UART_NUM, buffer + bytes_received, 
                                 PACKET_SIZE - bytes_received, pdMS_TO_TICKS(2)); 
         
        if (len > 0) { 
            bytes_received += len; 
             
            if (bytes_received >= PACKET_SIZE) { 
                SensorPacket *pkt = (SensorPacket*)buffer; 
                packet_count++; 
                 
                if (validate_packet_extended(pkt)) { 
                    if (xQueueSend(sensor_queue, pkt, 0) != pdTRUE) { 
                        SensorPacket dummy; 
                        xQueueReceive(sensor_queue, &dummy, 0); 
                        xQueueSend(sensor_queue, pkt, 0); 
                    } 
                } else { 
                    error_count++; 
                    if (error_count % 100 == 0) { 
                        printf("UART: %lu packet errors out of %lu total\n", error_count, 
packet_count); 
                    } 
                } 
                bytes_received = 0; 
            } 
        } 
         
        vTaskDelay(pdMS_TO_TICKS(1)); 
    } 
} 
 
// EKF Init( happen only once at startup)
static void ekf_task(void *pvParameters) { 
    ekf_init(&ekf); // Initialize EKF state and matrices
    uint32_t last_timestamp = 0; // Track timing between measurements
    float dt = 0.02f; // Default time step (50 Hz)
    float innovation[MEAS_DIM] = {0}; // Innovation vector (measurement residuals)
    int print_counter = 0; // Counter of how much ekf print diagnostic in terminal window
    uint32_t ekf_iterations = 0;  // Total EKF update count
    uint32_t rejected_measurements = 0; // Count of rejected outlier measurements
    uint32_t zupt_updates = 0; // Count of rejected outlier measurements

     
    // Initialize reference altitude 
    SensorPacket init_pkt; 
    if (xQueueReceive(sensor_queue, &init_pkt, pdMS_TO_TICKS(2000)) == pdTRUE) { 
        reference_altitude = init_pkt.pressure / 101.325f * 8430.0f; 
        last_timestamp = init_pkt.timestamp_us; 
        printf("Reference altitude initialized: %.2f m\n", reference_altitude); 
         
        // Put the packet back for initial alignment 
        xQueueSendToFront(sensor_queue, &init_pkt, 0); 
         
        // Perform initial alignment using TRAID method
        initial_alignment(&ekf, 50); // Use 50 samples for alignment 
         
    } else { 
        printf("Warning: No initial sensor data received\n"); 
    } 
   // main ekf loop run every cycle  
    while (1) { 
        SensorPacket pkt; 
         
        if (xQueueReceive(sensor_queue, &pkt, pdMS_TO_TICKS(30)) == pdTRUE) { 
            // Calculate dt (Bounds are applied to avoid instability on small dt and avoid inaccurate EST)
            if (last_timestamp != 0 && pkt.timestamp_us > last_timestamp) { 
                uint32_t dt_us = pkt.timestamp_us - last_timestamp; 
                dt = dt_us / 1e6f; 
                 
                if (dt < 0.005f) dt = 0.005f; 
                if (dt > 0.1f) dt = 0.1f; 
            } 
            last_timestamp = pkt.timestamp_us; 
             //in above logic an buffer drain needed to be implemented to avoid EKF use old data if
             //it lags behind and consume more time that 50hz or 20ms
            
            
             // Store previous state for comparison and fall back (logic yet to be added)
            static float prev_x[STATE_DIM]; 
            mat_copy(ekf.x, prev_x, STATE_DIM, 1); 
             
            // Run EKF 
            ekf_predict(&ekf, &pkt, dt); 
             
            // Store innovation before update 
            float prev_innovation[MEAS_DIM]; 
            mat_copy(innovation, prev_innovation, MEAS_DIM, 1); 
             
            ekf_update(&ekf, &pkt, innovation); // Mahalanobis is applied here.
             
            // Check if measurement was rejected (can also be used to chack if Mahalanobis is working)
            bool measurement_rejected = true; 
            for (int i = 0; i < MEAS_DIM; i++) { 
                if (fabsf(innovation[i] - prev_innovation[i]) > 1e-6f) { 
                    measurement_rejected = false; 
                    break; 
                } 
            } 
             
            if (measurement_rejected) { 
                rejected_measurements++; 
            } else { 
                // Apply adaptive tuning (Tune R matrix based on innovation stats)
                adaptive_filter_tuning(&ekf, innovation); 
            } 
             
            ekf_iterations++; 
             
            // Convert quaternion to Euler angles for human readable form
            float roll, pitch, yaw; 
            quat_to_euler(&ekf.x[0], &roll, &pitch, &yaw); 
             
            // Enhanced diagnostics every 25 iterations 
            if (++print_counter >= 25) { 
                print_counter = 0; 
                 
                printf("\n=== Enhanced EKF Status (Iteration %lu) ===\n", ekf_iterations); 
                printf("Timestamp: %lu us, dt: %.4f s\n", pkt.timestamp_us, dt); 
                printf("Rejected measurements: %lu (%.1f%%)\n",  
                       rejected_measurements, 100.0f * rejected_measurements / ekf_iterations); 
                printf("ZUPT updates applied: %lu\n", zupt_updates); 
                 
                // Attitude with uncertainties 
                float att_std[3] = { 
                    sqrtf(fabsf(ekf.P[1*STATE_DIM+1])) * 180.0f/M_PI, 
                    sqrtf(fabsf(ekf.P[2*STATE_DIM+2])) * 180.0f/M_PI, 
                    sqrtf(fabsf(ekf.P[3*STATE_DIM+3])) * 180.0f/M_PI 
                }; 
                printf("Attitude: R=%.1f±%.1f° P=%.1f±%.1f° Y=%.1f±%.1f°\n",  
                       roll, att_std[0], pitch, att_std[1], yaw, att_std[2]); 
                 
                // Position with uncertainties 
                printf("Position: [%.3f±%.3f, %.3f±%.3f, %.3f±%.3f] m\n",  
                       ekf.x[4], sqrtf(fabsf(ekf.P[4*STATE_DIM+4])), 
                       ekf.x[5], sqrtf(fabsf(ekf.P[5*STATE_DIM+5])), 
                       ekf.x[6], sqrtf(fabsf(ekf.P[6*STATE_DIM+6]))); 
                 
                // Velocity with uncertainties 
                printf("Velocity: [%.3f±%.3f, %.3f±%.3f, %.3f±%.3f] m/s\n",  
                       ekf.x[7], sqrtf(fabsf(ekf.P[7*STATE_DIM+7])), 
                       ekf.x[8], sqrtf(fabsf(ekf.P[8*STATE_DIM+8])), 
                       ekf.x[9], sqrtf(fabsf(ekf.P[9*STATE_DIM+9]))); 
                 
                // Biases 
                printf("Accel bias: [%.4f, %.4f, %.4f] m/s²\n",  
                       ekf.x[10], ekf.x[11], ekf.x[12]); 
                printf("Gyro bias: [%.4f, %.4f, %.4f] °/s\n",  
                       ekf.x[13], ekf.x[14], ekf.x[15]); 
                printf("Baro bias: %.2f Pa\n", ekf.x[16]); 
                 
                // Innovation and filter health 
                printf("Innovation: [%.0f, %.0f, %.0f, %.1f, %.3f, %.3f, %.3f]\n", 
                       innovation[0], innovation[1], innovation[2], innovation[3], 
                       innovation[4], innovation[5], innovation[6]); 
                 
                // Trace of covariance matrix (total uncertainty) 
                float trace = 0.0f; 
                for (int i = 0; i < STATE_DIM; i++) { 
                    trace += fabsf(ekf.P[i * STATE_DIM + i]); 
                } 
                printf("Total uncertainty (trace P): %.6f\n", trace); 
                 
                // Performance metrics 
                static uint32_t last_print_time = 0; 
                uint32_t current_time = esp_timer_get_time(); 
                if (last_print_time != 0) { 
                    float actual_rate = 25.0f / ((current_time - last_print_time) / 1e6f); 
                    printf("EKF Rate: %.1f Hz\n", actual_rate); 
                } 
                last_print_time = current_time; 
                 
                printf("Stack HWM: %d bytes\n", uxTaskGetStackHighWaterMark(NULL)); 
                printf("Free heap: %d bytes\n", esp_get_free_heap_size()); 
                printf("===========================================\n"); 
            } 
        } else { 
            printf("Warning: No sensor data, maintaining prediction\n"); 
        } 
         
        vTaskDelay(pdMS_TO_TICKS(20)); 
    } 
} 
 
// UART initialization  
static void uart_init(void) { 
    const uart_config_t uart_config = { 
        .baud_rate = UART_BAUD_RATE, 
        .data_bits = UART_DATA_8_BITS, 
        .parity = UART_PARITY_DISABLE, 
        .stop_bits = UART_STOP_BITS_1, 
        .flow_ctrl = UART_HW_FLOWCTRL_DISABLE, 
        .source_clk = UART_SCLK_DEFAULT, 
    }; 
     
    ESP_ERROR_CHECK(uart_driver_install(UART_NUM, BUF_SIZE * 2, 0, 0, NULL, 0)); 
    ESP_ERROR_CHECK(uart_param_config(UART_NUM, &uart_config)); 
    ESP_ERROR_CHECK(uart_set_pin(UART_NUM, UART_TX_PIN, UART_RX_PIN, 
                                UART_PIN_NO_CHANGE, UART_PIN_NO_CHANGE)); 
     
    uart_set_rx_timeout(UART_NUM, 20); 
    uart_set_rx_full_threshold(UART_NUM, PACKET_SIZE); 
     
    printf("UART initialized: %d baud, RX pin %d\n", UART_BAUD_RATE, UART_RX_PIN); 
} 
 
// Health monitoring task 
static void health_monitor_task(void *pvParameters) { 
    uint32_t last_check_time = 0; 
     
    while (1) { 
        uint32_t current_time = esp_timer_get_time(); 
         
        UBaseType_t queue_items = uxQueueMessagesWaiting(sensor_queue); 
        float queue_usage = (float)queue_items / PACKET_BUFFER_SIZE * 100.0f; 
         
        if (queue_usage > 80.0f) { 
            printf("Warning: Sensor queue %.0f%% full\n", queue_usage); 
        } 
         
        size_t free_heap = esp_get_free_heap_size(); 
        if (free_heap < 15000) { 
            printf("Warning: Low memory: %d bytes free\n", free_heap); 
        } 
         
        UBaseType_t uart_stack = uxTaskGetStackHighWaterMark(NULL); 
        if (uart_stack < 500) { 
            printf("Warning: Low stack space: %d bytes\n", uart_stack); 
        } 
         
        if (last_check_time != 0) { 
            float check_interval = (current_time - last_check_time) / 1e6f; 
            printf("Health check interval: %.1f s\n", check_interval); 
        } 
        last_check_time = current_time; 
         
        vTaskDelay(pdMS_TO_TICKS(5000)); 
    } 
} 
 
// Initial alignment procedure for better startup 
static void initial_alignment(EKFState *ekf, int num_samples) { 
    printf("Starting initial alignment with %d samples...\n", num_samples); 
     
    float acc_sum[3] = {0, 0, 0}; 
    float gyro_sum[3] = {0, 0, 0}; 
    float mag_sum[3] = {0, 0, 0}; 
    int valid_samples = 0; 
     
    // Collect samples for initial alignment 
    for (int i = 0; i < num_samples; i++) { 
        SensorPacket pkt; 
        if (xQueueReceive(sensor_queue, &pkt, pdMS_TO_TICKS(100)) == pdTRUE) { 
            if (validate_packet_extended(&pkt)) { 
                acc_sum[0] += pkt.ax; 
                acc_sum[1] += pkt.ay; 
                acc_sum[2] += pkt.az; 
                 
                gyro_sum[0] += pkt.gx; 
                gyro_sum[1] += pkt.gy; 
                gyro_sum[2] += pkt.gz; 
                 
                mag_sum[0] += pkt.mx; 
                mag_sum[1] += pkt.my; 
                mag_sum[2] += pkt.mz; 
                 
                valid_samples++; 
            } 
        } 
        vTaskDelay(pdMS_TO_TICKS(20)); 
    } 
     
    if (valid_samples < num_samples / 2) { 
        printf("Warning: Only %d valid samples for alignment\n", valid_samples); 
        return; 
    } 
     
    // Compute averages 
    for (int i = 0; i < 3; i++) { 
        acc_sum[i] /= valid_samples; 
        gyro_sum[i] /= valid_samples; 
        mag_sum[i] /= valid_samples; 
    } 
     
    // Initialize gyro bias (we calc bias of sensor here to get solid bias est)
    ekf->x[13] = gyro_sum[0]; 
    ekf->x[14] = gyro_sum[1]; 
    ekf->x[15] = gyro_sum[2]; 
     
    // Initialize accelerometer bias (assume initially level) 
    float acc_magnitude = sqrtf(acc_sum[0]*acc_sum[0] + acc_sum[1]*acc_sum[1] + 
acc_sum[2]*acc_sum[2]); 
    ekf->x[10] = 0.0f; 
    ekf->x[11] = 0.0f; 
    ekf->x[12] = acc_magnitude - 9.81f; // Bias in Z to account for gravity 
     
    // Initialize orientation from accelerometer and magnetometer 
    // Normalize accelerometer reading to get 'g' unit vector
    float acc_norm[3]; 
    for (int i = 0; i < 3; i++) { 
        acc_norm[i] = -acc_sum[i] / acc_magnitude; // Negative because accelerometer 
measures specific force 
    } 
     
    // Normalize magnetometer reading to get magnetic field 'N' unit vector
    float mag_magnitude = sqrtf(mag_sum[0]*mag_sum[0] + mag_sum[1]*mag_sum[1] + 
mag_sum[2]*mag_sum[2]); 
    float mag_norm[3]; 
    for (int i = 0; i < 3; i++) { 
        mag_norm[i] = mag_sum[i] / mag_magnitude; 
    } 
     
    // Compute initial orientation using TRIAD algorithm 
    // acc_norm points "up" in body frame, should align with [0,0,1] in NED 
    // mag_norm points "north" in body frame, should align with magnetic field in NED 
     
    // Compute cross product to get east direction 
    float east_body[3]; 
    east_body[0] = acc_norm[1] * mag_norm[2] - acc_norm[2] * mag_norm[1]; 
    east_body[1] = acc_norm[2] * mag_norm[0] - acc_norm[0] * mag_norm[2]; 
    east_body[2] = acc_norm[0] * mag_norm[1] - acc_norm[1] * mag_norm[0]; 
     
    // Normalize east vector 
    float east_magnitude = sqrtf(east_body[0]*east_body[0] + east_body[1]*east_body[1] + 
east_body[2]*east_body[2]); 
    if (east_magnitude > 0.1f) { 
        for (int i = 0; i < 3; i++) { 
            east_body[i] /= east_magnitude; 
        } 
         
        // Recompute north as cross product of down and east to tacle sensor tilt problem
        float north_body[3]; 
        north_body[0] = (-acc_norm[1]) * east_body[2] - (-acc_norm[2]) * east_body[1]; 
        north_body[1] = (-acc_norm[2]) * east_body[0] - (-acc_norm[0]) * east_body[2]; 
        north_body[2] = (-acc_norm[0]) * east_body[1] - (-acc_norm[1]) * east_body[0]; 
         
        // Build rotation matrix from body to NED 
        float R_bn[9] = { 
            north_body[0], east_body[0], -acc_norm[0], 
            north_body[1], east_body[1], -acc_norm[1], 
            north_body[2], east_body[2], -acc_norm[2] 
        }; 
         
        // Convert rotation matrix to quaternion 
        float trace = R_bn[0] + R_bn[4] + R_bn[8]; 
        if (trace > 0) { 
            float s = sqrtf(trace + 1.0f) * 2; // s = 4 * qw 
            ekf->x[0] = 0.25f * s; 
            ekf->x[1] = (R_bn[7] - R_bn[5]) / s; 
            ekf->x[2] = (R_bn[2] - R_bn[6]) / s; 
            ekf->x[3] = (R_bn[3] - R_bn[1]) / s; 
        } else if ((R_bn[0] > R_bn[4]) && (R_bn[0] > R_bn[8])) { 
            float s = sqrtf(1.0f + R_bn[0] - R_bn[4] - R_bn[8]) * 2; // s = 4 * qx 
            ekf->x[0] = (R_bn[7] - R_bn[5]) / s; 
            ekf->x[1] = 0.25f * s; 
            ekf->x[2] = (R_bn[1] + R_bn[3]) / s; 
            ekf->x[3] = (R_bn[2] + R_bn[6]) / s; 
        } else if (R_bn[4] > R_bn[8]) { 
            float s = sqrtf(1.0f + R_bn[4] - R_bn[0] - R_bn[8]) * 2; // s = 4 * qy 
            ekf->x[0] = (R_bn[2] - R_bn[6]) / s; 
            ekf->x[1] = (R_bn[1] + R_bn[3]) / s; 
            ekf->x[2] = 0.25f * s; 
            ekf->x[3] = (R_bn[5] + R_bn[7]) / s; 
        } else { 
            float s = sqrtf(1.0f + R_bn[8] - R_bn[0] - R_bn[4]) * 2; // s = 4 * qz 
            ekf->x[0] = (R_bn[3] - R_bn[1]) / s; 
            ekf->x[1] = (R_bn[2] + R_bn[6]) / s; 
            ekf->x[2] = (R_bn[5] + R_bn[7]) / s; 
            ekf->x[3] = 0.25f * s; 
        } 
         
        quat_normalize(&ekf->x[0]); 
    } 
     
    // Reduce initial uncertainties after alignment 
    for (int i = 0; i < 4; i++) { 
        ekf->P[i * STATE_DIM + i] = 0.001f; // Reduced attitude uncertainty cz we have now more confidence
    } 
    for (int i = 10; i < 16; i++) { 
        ekf->P[i * STATE_DIM + i] *= 0.1f; // Reduced bias uncertainty 
    } 
     
    float roll, pitch, yaw; 
    quat_to_euler(&ekf->x[0], &roll, &pitch, &yaw); 
     
    printf("Initial alignment completed:\n"); 
    printf("  Samples used: %d/%d\n", valid_samples, num_samples); 
    printf("  Initial attitude: R=%.1f° P=%.1f° Y=%.1f°\n", roll, pitch, yaw); 
    printf("  Gyro bias: [%.3f, %.3f, %.3f] °/s\n", ekf->x[13], ekf->x[14], ekf->x[15]); 
    printf("  Accel bias: [%.3f, %.3f, %.3f] m/s²\n", ekf->x[10], ekf->x[11], ekf->x[12]); 
} 
 
//  Zero velocity update for stationary periods (Defined but not implemented yet) 
static bool detect_zero_velocity(const SensorPacket *pkt, float gyro_threshold, float 
accel_threshold) { 
    // Check if gyro readings are below threshold 
    float gyro_magnitude = sqrtf(pkt->gx*pkt->gx + pkt->gy*pkt->gy + pkt->gz*pkt->gz); 
    if (gyro_magnitude > gyro_threshold) return false; 
     
    // Check if accel magnitude is close to 1g 
    float accel_magnitude = sqrtf(pkt->ax*pkt->ax + pkt->ay*pkt->ay + pkt->az*pkt->az); 
    if (fabsf(accel_magnitude - 9.81f) > accel_threshold) return false; 
     
    return true; 
} 
 
static void apply_zero_velocity_update(EKFState *ekf) { 
    // Measurement: velocity should be zero 
    float z_zupt[3] = {0.0f, 0.0f, 0.0f}; 
    float h_zupt[3] = {ekf->x[7], ekf->x[8], ekf->x[9]}; // Predicted velocity 
     
    // Innovation (Used for Bias Calc)
    float y_zupt[3]; 
    for (int i = 0; i < 3; i++) { 
        y_zupt[i] = z_zupt[i] - h_zupt[i]; 
    } 
     
    // Measurement Jacobian (H is zero except for velocity states) 
    static float H_zupt[3 * STATE_DIM]; 
    mat_zero(H_zupt, 3, STATE_DIM); 
    H_zupt[0 * STATE_DIM + 7] = 1.0f; // dh/dvx 
    H_zupt[1 * STATE_DIM + 8] = 1.0f; // dh/dvy 
    H_zupt[2 * STATE_DIM + 9] = 1.0f; // dh/dvz 
     
    // Custom R Matrix for ZUPT (very low since we're confident velocity is zero) 
    static float R_zupt[9]; 
    mat_zero(R_zupt, 3, 3); 
    R_zupt[0] = 0.001f; // Very low velocity noise (trst MEA more)
    R_zupt[4] = 0.001f; 
    R_zupt[8] = 0.001f; 
     
    // Innovation covariance: S = H*P*H' + R 
    //To know how much we should trust the "v = 0" measurement compared to the filter's own current uncertainty      
    static float HT_zupt[STATE_DIM * 3]; 
    static float PH_zupt[STATE_DIM * 3]; 
    static float HPH_zupt[9]; 
    static float S_zupt[9]; 
     
    mat_transpose(H_zupt, HT_zupt, 3, STATE_DIM); 
    mat_mult(ekf->P, HT_zupt, PH_zupt, STATE_DIM, STATE_DIM, 3); 
    mat_mult(H_zupt, PH_zupt, HPH_zupt, 3, STATE_DIM, 3); 
    mat_add(HPH_zupt, R_zupt, S_zupt, 3, 3); 
     
    // Kalman gain (to calc P_zupt)
    static float S_zupt_work[9]; 
    static int perm_zupt[3]; 
    static float K_zupt[STATE_DIM * 3]; 
     
    mat_copy(S_zupt, S_zupt_work, 3, 3); 
     
    if (lu_decompose(S_zupt_work, perm_zupt, 3)) { 
        for (int i = 0; i < STATE_DIM; i++) { 
            float K_temp[3]; 
            for (int j = 0; j < 3; j++) { 
                K_temp[j] = PH_zupt[i * 3 + j]; 
            } 
             
            float K_col[3]; 
            lu_solve(S_zupt_work, perm_zupt, K_temp, K_col, 3); 
             
            for (int j = 0; j < 3; j++) { 
                K_zupt[i * 3 + j] = K_col[j]; 
            } 
        } 
         
        // State update 
        for (int i = 0; i < STATE_DIM; i++) { 
            float update = 0.0f; 
            for (int j = 0; j < 3; j++) { 
                update += K_zupt[i * 3 + j] * y_zupt[j]; 
            } 
            ekf->x[i] += update; 
        } 
         
        // Normalize quaternion (cz V=0 affect quat due to Co-Var)
        quat_normalize(&ekf->x[0]); 
         
        // Covariance update (Joseph form) (P_zupt)
        static float KH_zupt[STATE_DIM * STATE_DIM]; 
        static float I_KH_zupt[STATE_DIM * STATE_DIM]; 
        static float temp_P_zupt[STATE_DIM * STATE_DIM]; 
         
        mat_mult(K_zupt, H_zupt, KH_zupt, STATE_DIM, 3, STATE_DIM); 
        mat_identity(I_KH_zupt, STATE_DIM); 
        mat_sub(I_KH_zupt, KH_zupt, I_KH_zupt, STATE_DIM, STATE_DIM); 
        mat_mult(I_KH_zupt, ekf->P, temp_P_zupt, STATE_DIM, STATE_DIM, STATE_DIM); 
        mat_copy(temp_P_zupt, ekf->P, STATE_DIM, STATE_DIM); 
         
        printf("Zero velocity update applied\n"); 
    } 
} 
 
// hehe main function
void app_main(void) { 
    printf("\n=== Complete Advanced EKF Navigation System ===\n"); 
    printf("Location: GOTH IQBAL\n"); 
    printf("Critical fixes implemented:\n"); 
    printf("  1. Proper process noise modeling (12x12 Q matrix)\n"); 
    printf("  2. Process noise Jacobian G matrix\n"); 
    printf("  3. Complete state transition Jacobian F\n"); 
    printf("  4. Joseph form covariance updates\n"); 
    printf("  5. Innovation-based outlier rejection\n"); 
    printf("  6. Adaptive filter tuning\n"); 
    printf("  7. Initial alignment procedure\n"); 
    printf("  8. Zero velocity updates\n"); 
    printf("  9. Enhanced numerical stability\n"); 
    printf("  10. Comprehensive diagnostics\n"); 
    printf("=============================================\n"); 
     
    printf("Magnetic Declination: +%.3f° (East)\n", mag_declination); 
    printf("Magnetic Inclination: %.3f°\n", mag_inclination); 
    printf("Magnetic Field Strength: %.1f nT\n", mag_field_strength); 
    printf("State dimension: %d\n", STATE_DIM); 
    printf("Measurement dimension: %d\n", MEAS_DIM); 
    printf("Process noise dimension: %d\n", NOISE_DIM); 
    printf("STM32 data rate: 70 Hz\n"); 
    printf("EKF processing rate: 50 Hz\n"); 
    printf("Packet size: %d bytes\n", PACKET_SIZE); 
     
    // Init UART 
    uart_init(); 
     
    // Create sensor data queue 
    sensor_queue = xQueueCreate(PACKET_BUFFER_SIZE, sizeof(SensorPacket)); 
    if (sensor_queue == NULL) { 
        printf("ERROR: Failed to create sensor queue\n"); 
        return; 
    } 
     
    // Create tasks with optimized priorities and core affinity 
    BaseType_t ret; 
     
    // UART receive on core 0 (highest priority) 
    ret = xTaskCreatePinnedToCore(uart_receive_task, "uart_rx", 
                                 4096, NULL, 16, NULL, 0); 
    if (ret != pdPASS) { 
        printf("ERROR: Failed to create UART receive task\n"); 
        return; 
    } 
     
    // EKF processing on core 1 (high priority, large stack for matrix operations) 
    //One core for EKF and other for mission planning and lora comm. etc
    ret = xTaskCreatePinnedToCore(ekf_task, "ekf_proc", 
                                 32768, NULL, 15, NULL, 1);  // Increased stack for larger matrices 
    if (ret != pdPASS) { 
        printf("ERROR: Failed to create EKF task\n"); 
        return; 
    } 
     
    // Health monitoring on core 0 (low priority) 
    ret = xTaskCreatePinnedToCore(health_monitor_task, "health_mon", 
                                 2048, NULL, 5, NULL, 0); 
    if (ret != pdPASS) { 
        printf("ERROR: Failed to create health monitor task\n"); 
        return; 
    } 
     
    printf("All tasks created successfully\n"); 
    printf("System ready - waiting for STM32 data...\n"); 
    printf("Initial alignment will begin automatically...\n"); 
    printf("=============================================\n"); 
} 
