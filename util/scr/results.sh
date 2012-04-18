#!/bin/bash
fracen=$(grep "Fractional energy" information.out | grep "integrator" | tail -n 1 | cut -f 2 -d ':')
fracmom=$(grep "Fractional angular" information.out | tail -n 1 | cut -f 2 -d ':')
{
    echo "ChangeEnergy=" $fracen
    echo "ChangeMomentum=" $fracmom
} > sci2web/results.info
