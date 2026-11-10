#ifndef AMOSTRAGEM_H
#define AMOSTRAGEM_H
#include <cstdint>

namespace amostragem {
// J e P extraidas de benchmarking.cpp, suplemento de Ahmed (EGSR 2024).
inline std::uint32_t J(std::uint32_t i) {
    i = (i << 16) | (i >> 16);
    i = ((i & 0x00FF00FFu) << 8) | ((i & 0xFF00FF00u) >> 8);
    i = ((i & 0x0F0F0F0Fu) << 4) | ((i & 0xF0F0F0F0u) >> 4);
    i = ((i & 0x33333333u) << 2) | ((i & 0xCCCCCCCCu) >> 2);
    i = ((i & 0x55555555u) << 1) | ((i & 0xAAAAAAAAu) >> 1);
    return i;
}
inline std::uint32_t P(std::uint32_t v) {
    v ^= v << 16;
    v ^= (v & 0x00FF00FFu) << 8;
    v ^= (v & 0x0F0F0F0Fu) << 4;
    v ^= (v & 0x33333333u) << 2;
    v ^= (v & 0x55555555u) << 1;
    return v;
}
// Mistura deterministica para sementes. Nao faz parte do algoritmo de Ahmed.
inline std::uint32_t mistura(std::uint32_t v) {
    v ^= v >> 16; v *= 0x7feb352du;
    v ^= v >> 15; v *= 0x846ca68bu;
    return v ^ (v >> 16);
}
inline void sobol(std::uint32_t indice, std::uint32_t chave,
                  double &u, double &v) {
    const std::uint32_t x = J(indice);
    // Um deslocamento digital fixo por subpixel evita repetir o mesmo
    // padrao em toda a imagem. Nao e Owen scrambling.
    u = double(x ^ mistura(chave ^ 0x68bc21ebu)) / 4294967296.0;
    v = double(P(x) ^ mistura(chave ^ 0x02e5be93u)) / 4294967296.0;
}
inline void estado(std::uint32_t chave, unsigned short xi[3]) {
    const std::uint32_t a = mistura(chave), b = mistura(chave ^ 0xa511e9b3u);
    xi[0] = static_cast<unsigned short>(a);
    xi[1] = static_cast<unsigned short>(a >> 16);
    xi[2] = static_cast<unsigned short>(b);
}
}
#endif
