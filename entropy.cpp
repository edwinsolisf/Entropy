#include <iostream>
#include <cmath>
#include <string>
#include <cstdio>

double temperature(uint64_t atoms, uint64_t quanta, double stiffness, double mass);
double specific_heat(uint64_t atoms, uint64_t quanta, double stiffness, double mass);

/**
 * @args numberOfAtoms stiffness(N/m) massOfAtom(kg)
 * Produces a text file with data order in Atoms-Quanta-Temperature-Cv
 * This file can be imported into excel
 */
int main(int argc, char** argv)
{
    int atoms = 1;
    double stiffness = 1;
    double mass = 1;
    if (argc > 3)
    {
        atoms = std::stoi(argv[1]);
        stiffness = std::stod(argv[2]);
        mass = std::stod(argv[3]);
    }

    FILE* file = fopen("data.txt", "w");

    fprintf(file, "Atoms\tQuanta\tTemperature\tCv\n");


    for (int i = 1; i < 100; ++i)
        fprintf(file, "%d\t%d\t%E\t%E\n", atoms, i, temperature(atoms, i, stiffness, mass), specific_heat(atoms, i, stiffness, mass));
    for (int i = 100; i < 10000; i += 100)
        fprintf(file, "%d\t%d\t%E\t%E\n", atoms, i, temperature(atoms, i, stiffness, mass), specific_heat(atoms, i, stiffness, mass));
    for (int i = 10000; i < 1000000; i += 10000)
        fprintf(file, "%d\t%d\t%E\t%E\n", atoms, i, temperature(atoms, i, stiffness, mass), specific_heat(atoms, i, stiffness, mass));
    for (int i = 1000000; i < 10000000; i += 100000)
        fprintf(file, "%d\t%d\t%E\t%E\n", atoms, i, temperature(atoms, i, stiffness, mass), specific_heat(atoms, i, stiffness, mass));

    fclose(file);

    return 0;
}

double specific_heat(uint64_t atoms, uint64_t quanta, double stiffness, double mass)
{
    constexpr double R = 8.31446261815324;
    uint64_t oscillators = atoms * 3;

    double diffs32 = log((double)(oscillators + quanta)/(double)(quanta + 1));
    double diffs21 = log((double)(oscillators + quanta - 1)/(double)(quanta));

    return R * (diffs32 * diffs21 / (diffs21 - diffs32)) / atoms;
}

double temperature(uint64_t atoms, uint64_t quanta, double stiffness, double mass)
{
    constexpr double kb = 1.380649E-23;
    constexpr double hbar = 6.62607015E-34 / (2.0 * 3.14159265358979323846);

    uint64_t oscillators = atoms * 3;
    double diffs32 = log((double)(oscillators + quanta)/(double)(quanta + 1));
    double diffs21 = log((double)(oscillators + quanta - 1)/(double)(quanta));
    double diffs31 = diffs32 + diffs21;
    double diffomega = sqrt(4.0 * stiffness / mass);

    double change = diffs31 / 2;

    return (hbar * diffomega) / (kb * change);
}
