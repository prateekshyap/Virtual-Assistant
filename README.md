# Virtual Assistant

This is a command based version of Yugi. For the GUI version please look at [Yugi](https://github.com/therohanjaiswal/Yugi)

# Steps to execute

1. Open [FinalProject.sln](https://github.com/prateekshyap/Tutorial-for-kids/blob/main/FinalProject.sln) file in Visual Studio 2010.
2. Press F5 or ctrl+F5 for building and running the solution.
3. You can start with the existing model.
4. If you want to go back to the new-born version, store __0__ in [info.txt](https://github.com/prateekshyap/Tutorial-for-kids/blob/main/FinalProject/data/info.txt). It will delete all extra trainings and will keep only the default words which were added in the beginning of this project.
5. Then you can try live testing. Speak the word within 3 seconds. The resultant word will be displayed. (Browsing feature is not added here. It is present in the GUI version.)
6. For live training, first write the word you want to train. If you're adding a new word, you'll have to record 10 times continuously.
7. Following the switch case and instructions will be enough.

## Suggestions or improvements?

Fork this repository, make required changes and create a pull request.


----------------------------------------------------------------------------------------------------------
# Past information (not for viewers)

## command to record-
```Recording_Module.exe 3 WAVFiles/car/0.wav FinalProject/HMM/car/0.txt```

## Recordings-

### Recorded - 

null<br>
aeroplane<br>
apple<br>
autorickshaw<br>
bicycle<br>
bike<br>
bus<br>
car<br>

### Pending recordings -

cat<br>
dog<br>
scooter<br>
tandem<br>
train<br>
tram<br>
truck<br>

Please record the __Recorded__ words first, in the above sequence, 0.txt to 9.txt. 10.txt to 19.txt already added.

## Different header files to be included-

1. Vowel Recognition Assignment - RJ (implemented)
2. LBG - PP (implemented)
3. Forward and Backward Process - RJ (implemented)
4. Viterbi Algorithm - PP (implemented)
5. Baum Welch Method - PP (implemented)

## Current status

~~This code is now able to take input for A, B, O and Pi and generate A', B' and Pi' for the given O.~~

~~This code is now able to take input for A, B, O and Pi and converge the model to an optimal one for each observation sequence.~~

~~This code is now able to train a model for the given words and recognize a single word spoken by an user.~~

Command based application for Yugi. For the GUI version please look at [Yugi](https://github.com/therohanjaiswal/Yugi)

## __Please add image display feature or youtube link feature for any given word.__
The code snippet should be able to use a ``` char * ``` and display some image or video related to it.<br>
Run the code in a different cpp file and paste the code in ```FinalProject.cpp```.<br>
Do not modify existing code.
