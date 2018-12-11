# MPC-with-PHE

Proof of concept code for implementing a Model Predictive Control problem on encrypted data in Python. More details can be found in the paper "Cloud-based MPC with Encrypted Data" https://arxiv.org/pdf/1803.09891.pdf.

This code is intended to give an idea on the execution speed and memory, and not for application that require highly optimized code.

Both the code for the Client-Server architecture and Two-Server architecture are added.

The file runCS.py and runSS.py are the files that run a thread for each of the two computing parties. 

The matrices and vectors for the Model Predictive Control formulation as in the paper have to be written in files in the folder Data, with the name having name-of-matrix_number-of-states_number-of-inputs_horizon.txt. 

The Paillier (and DGK for the two server case) keys should be written in files in the folder Keys. If this is not desired, the code has to be manually modified to generate them. The code will be updated soon to perform the check automatically.

The random coins for encryption and blinding can be pre-generated in order to save time in the online execution, and stored in a folder Randomness. At the moment, the code checks if files with random numbers of appropiate size exists, otherwise, it generates them from scratch.

Comments in the code will be added in the near future. A tidier version will also be added in the near future.
