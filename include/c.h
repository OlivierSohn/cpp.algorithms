#ifndef IMJ_AUDIO_C_H
#define IMJ_AUDIO_C_H

/********************************************/
// These structs are reflected in Haskell in .hsc files.
// The types of the members must match (in size)
//   their counterparts in Haskell data records.
//   For example we associate 'long' with Int32, and float with 'Float'.
//   We don't use 'int' because its size may vary depending on the platform.
//
// Hence, when modifying a type here, please make sure to
// modify also the type in the corresponding .hsc file.
/*********************************************/

typedef struct {

  // Expressed in 'radians / pi',
  // i.e [-1,1] covers the whole trigonometric circle
  float phase;

  // Between 0 (included) and 1 (included).
  float volume;
} harmonicProperties_t;

typedef struct {
  long nChannels;
  long sampleSize;
  long nFrames;
  long sampleRate;
  float lengthInSeconds;
} spaceResponse_t;

#endif
