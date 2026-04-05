# SVD & Solveurs Numériques

## Description
Implémentation de l'algorithme de Décomposition en Valeurs Singulières (SVD), pilier du pipeline AR pour la calibration et l'estimation de pose.

Ce projet fait partie du cursus de Réalité Augmentée From Scratch (Semestre 1). L'implémentation est réalisée entièrement en C++17 sans bibliothèque tierce (Zero STL), en suivant les directives de développement moteur Jenga/Nkentseu.

## Fonctionnalités implémentées
- Algorithme SVD Jacobi one-sided pour matrices rectangulaires.
- Calcul des valeurs singulières.
- Résolution de systèmes Ax=0.
- Calcul de la pseudo-inverse matricielle.

## Installation et Compilation
Le projet utilise le système de build **Jenga**.

1. Assurez-vous d'avoir Jenga installé sur votre système.
2. Clonez le dépôt :
   ```bash
   git clone git@github.com:neussi/RA_TP5_SVD_et_Solveurs.git
   cd RA_TP5_SVD_et_Solveurs
   ```
3. Compilez le projet :
   ```bash
   jenga build
   ```
4. Exécutez le programme :
   ```bash
   jenga run TP5
   ```

## Résultats
Le programme génère des sorties dans la console et, le cas échéant, des fichiers images (.ppm) illustrant les concepts mathématiques et graphiques abordés.

## Auteur
**NEUSSI NJIETCHEU PATRICE EUGNE**
Matricule : 24P820
