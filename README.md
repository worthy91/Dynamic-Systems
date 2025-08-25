# Dynamic-Systems
# 17-State Extended Kalman Filter Navigation System

A real-time navigation system implementing a 17-state Extended Kalman Filter for autonomous UAV applications, processing IMU and magnetometer data at 50Hz on ESP32 with STM32 sensor acquisition.

## Overview

This project implements a comprehensive navigation solution using sensor fusion techniques commonly found in aerospace and robotics applications. The system processes inertial measurement unit (IMU) data, magnetometer readings, and barometric pressure to estimate vehicle attitude, position, velocity, and sensor biases in real-time.

## System Architecture

- **Primary Processor**: ESP32 (EKF processing at 50Hz)
- **Sensor Interface**: BluePill STM32 (sensor acquisition at 70Hz)
- **Communication**: UART data streaming
- **Target Application**: Autonomous UAV navigation

## State Vector (17 States)

1. **Quaternion Attitude** (4 states): `[qw, qx, qy, qz]`
2. **3D Position** (3 states): `[px, py, pz]` in NED frame
3. **3D Velocity** (3 states): `[vx, vy, vz]` in NED frame
4. **Accelerometer Biases** (3 states): `[bax, bay, baz]`
5. **Gyroscope Biases** (3 states): `[bgx, bgy, bgz]`
6. **Barometer Bias** (1 state): `[bp]`

## Sensor Suite

- **9-DOF IMU**: 3-axis accelerometer, gyroscope, magnetometer
- **Barometric Pressure Sensor**: Altitude estimation
- **Magnetic Field Modeling**: Location-specific declination and inclination

## Key Features

### ✅ Implemented
- Complete Jacobian matrices (F, G, H) for proper EKF implementation
- Joseph-form covariance updates for numerical stability
- Mahalanobis distance outlier rejection (Chi-squared test)
- TRIAD initial alignment algorithm for startup
- Adaptive filter tuning based on innovation statistics
- 4th-order Runge-Kutta quaternion integration
- Custom matrix operations (no external dependencies)
- Real-time diagnostics and health monitoring
- Magnetic field modeling for specific geographic location

### ⚠️ Current Limitations
- Fixed Q matrix (needs adaptive process noise)
- Raw mathematical implementations instead of optimized libraries
- ZUPT (Zero Velocity Updates) coded but not implemented
- Batch Mahalanobis updates instead of sequential processing
- H matrix limited to attitude-based innovations
- No GPS integration
- Missing buffer drain mechanism for EKF lag
- No fallback recovery for system failures
- Requires level and stationary startup
- No temperature compensation for gyro drift

## Hardware Requirements

### ESP32 Module
- Minimum 4MB flash
- Dual-core processor recommended
- UART1 for STM32 communication

### STM32 BluePill
- STM32F103C8T6 or similar
- Connected to IMU sensors via I2C/SPI
- UART output to ESP32

### Sensors
- 9-DOF IMU (e.g., MPU9250, ICM20948)
- Barometric pressure sensor (e.g., BMP280, MS5611)

## Installation

### Prerequisites
- ESP-IDF framework
- STM32 development environment
- Connected sensor hardware

### Building
```bash
# Clone repository
git clone [your-repo-url]
cd 17-state-ekf-navigation

# Configure ESP-IDF
idf.py menuconfig

# Build and flash
idf.py build
idf.py -p [PORT] flash monitor
```

## Configuration

### Magnetic Field Parameters
Update the magnetic field parameters for your location in the code:
```c
static float mag_declination = 1.317f;    // Magnetic declination in degrees
static float mag_inclination = 39.417f;   // Magnetic inclination in degrees  
static float mag_field_strength = 45681.5f; // Total field strength in nT
```

### UART Configuration
```c
#define UART_BAUD_RATE 115200
#define UART_RX_PIN 18
```

## Usage

1. **Power up the system** with sensors level and stationary
2. **Wait for initial alignment** (automatic TRIAD algorithm)
3. **Monitor real-time output** showing:
   - Attitude (roll, pitch, yaw) with uncertainties
   - Position and velocity estimates
   - Sensor bias estimates
   - Filter health diagnostics

### Sample Output
```
=== Enhanced EKF Status (Iteration 1250) ===
Timestamp: 25000000 us, dt: 0.0200 s
Attitude: R=-1.2±0.5° P=2.1±0.4° Y=45.3±0.8°
Position: [0.123±0.056, -0.045±0.034, 1.234±0.123] m
Velocity: [0.001±0.012, 0.003±0.015, -0.002±0.018] m/s
```

## Mathematical Foundation

The system implements a full nonlinear EKF with:
- **Prediction Step**: State propagation using IMU data
- **Update Step**: Correction using magnetometer and barometer
- **Numerical Stability**: Joseph-form covariance updates
- **Outlier Rejection**: Mahalanobis distance gating

### Process Model
```
x(k+1) = f(x(k), u(k)) + w(k)
```

### Measurement Model
```
z(k) = h(x(k)) + v(k)
```

Where measurements include magnetometer readings, barometric altitude, and accelerometer data.

## Performance

- **Processing Rate**: 50Hz EKF updates
- **Data Rate**: 70Hz sensor acquisition
- **Memory Usage**: ~32KB stack for matrix operations
- **Latency**: <20ms processing time per update

## Development Status

This is a **learning project** and **not production-ready**. The implementation prioritizes educational value and understanding of navigation algorithms over optimization and robustness.

### Future Improvements
- [ ] Adaptive process noise covariance
- [ ] GPS integration
- [ ] Sequential measurement processing
- [ ] Temperature compensation
- [ ] Robust fallback mechanisms
- [ ] Flight testing and validation

## Contributing

Contributions are welcome, especially:
- Bug fixes and numerical stability improvements
- Additional sensor integrations
- Performance optimizations
- Documentation improvements

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Mathematical foundations based on "Aided Navigation" by Jay A. Farrell
- TRIAD algorithm implementation from spacecraft attitude determination literature
- Joseph-form covariance update for numerical stability

## Contact

- LinkedIn: [linkedin.com/in/ali-iqbal-2004-ev](https://linkedin.com/in/ali-iqbal-2004-ev)
- Email: aaauif1234@gmail.com

---

**⚠️ Disclaimer**: This is an educational project. Use in safety-critical applications requires extensive testing and validation.
