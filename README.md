# stp_lib
stopping power libary

stopping power or dE/dx is to describe how a particle loses its energy when it is traveling in a medium. For example, it can be used to estimate how much energy a particle will deposit when it passes through your detector. In nuclear physics, the unit of dE/dx is MeV/(mg/cm2), since the beam normally in unit of MeV, and the target thickness is in unit of mg/cm2. <br> <br>

There are several programs can simulate dE/dx for all kinds of particles in a medium. For example, SRIM and LISE++. However, sometimes, I need to generate dE/dx data for internal use for my program, and reading stopping power from external files is not genuinely convenient. So I search the library that can help me generate dE/dx file, but I don't have any good luck to find it.  <br><br>

I come across that VIKAR program ( Virtual Instrumentation for Kinematics and Reactions, by Dr. S.D. Pain ) has its dE/dx generator, and I have the source code. Not only I would like to study how VIKAR's dE/dx data compares to those from LISE++, but also I would like to extract VIKAR's dE/dx generator. The VIKAR is written in Fortran, and some parts of code is in Fortran 77. The main part for dE/dx calculation is the subroutine ncdedx. I rewrite and simplify the ncdedx subroutine, and translate the code into a C++ class. It provides easy-to-use interface methods to get the dE/dx and stopping range data. 

++++++++++++++++++++++++++++++++++++++++++++++++++++++

please see <a href="http://peiluan-tai.com/programs/stp_lib.html"> here </a> for how to use it.
