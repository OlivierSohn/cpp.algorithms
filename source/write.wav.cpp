

using namespace imajuscule::audio;

void WAVWriter::DoUpdateFileHeader()
{
    header.subchunk2_size = n_sample_bytes_written;

    WriteData(&header.chunk_id[0], 4, 1);
    WriteData(&header.chunk_size, sizeof(int32_t), 1);
    WriteData(&header.format[0], 4, 1);
    WriteData(&header.subchunk1_id[0], 4, 1);
    WriteData(&header.subchunk1_size, sizeof(int32_t), 1);
    WriteData(&header.audio_format, sizeof(int16_t), 1);
    WriteData(&header.num_channels, sizeof(int16_t), 1);
    WriteData(&header.sample_rate, sizeof(int32_t), 1);
    WriteData(&header.byte_rate, sizeof(int32_t), 1);
    WriteData(&header.block_align, sizeof(int16_t), 1);
    WriteData(&header.bits_per_sample, sizeof(int16_t), 1);
    WriteData(&header.subchunk2_id[0], 4, 1);
    WriteData(&header.subchunk2_size, sizeof(int32_t), 1);
}
