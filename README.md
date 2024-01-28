# CS922 lab files

[![Try Build](https://github.com/Xanonymous-GitHub/cs922-lab/actions/workflows/build.yml/badge.svg)](https://github.com/Xanonymous-GitHub/cs922-lab/actions/workflows/build.yml)

## Lab 1

In this lab, we are going to take a brief look at C and C++. The two languages are quite similar in terms of syntax and
are often grouped together when discussing various programming languages. We’ll go through the basics and defining
characteristics of both languages, discuss the tools required to compile and run source code, we’ll run through a few
exercises, and at the end there are a few problems for you to test out your knowledge.
If you have done C and C++ before and know the basics including pointers memory allocation and classes, then you can
skip.

## Lab 2

OpenMP (the MP stands for Multiprocessing) is an API and runtime which enables the programming of multiple processing
cores with shared memory.
The API is a collection of functions and pragmas; the former allows the querying of information such as the number of
active threads,
and the latter allows the definition of parallel regions. In these parallel regions, threads are spawned by the OpenMP
runtime and the code within is executed on multiple cores.
These threads are then able to communicate with each other using shared memory.
This lab will give a basic introduction on how to use the OpenMP API to exploit parallelism in applications. Whilst
doing this lab,
it is recommended that you have the OpenMP Cheat Sheet.
The work in this lab requires the knowledge gained from Lab 1. As such, if you have not read/completed Sections 2 and 3
in Lab 1, then it is highly recommended you do this!

## DEQN (Coursework 1, lab3)

We are going to look at deqn, a simple simulation code that we will be using for Coursework 1.

The deqn code solves the diffusion equation in two-dimensions.

The diffusion equation is a partial differential equation that describes the transfer of heat through a material.

The deqn program is designed to solve the diffusion equation in two dimensions.
It is written completely in C++, and the program is encapsulated into multiple objects.
In this section we will walk through the source code of deqn, so that you can see how it all fits together.
