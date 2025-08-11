# Ubicacion Archivos

todos los archivos, fuentes y testbench deben estar en la
raiz de la carpeta del componentes de HLS, tomar como
ejemplo la estructura de Custom_Thomas_Solver_v1. 

esto facilita la declaraciones en hls_config.cfg

# Librería complex

- Usar siempre `std::complex` en lugar de `hls::complex`.
- `std::complex` es parte del estándar de C++ y tiene mejor compatibilidad; `hls::complex` queda solo por soporte heredado según AMD Xilinx Vitis HLS 2025.1.

