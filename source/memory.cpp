
#  ifdef IMJ_LOG_MEMORY

void* operator new  ( std::size_t count ) {
  auto ptr = std::malloc(count);

  imajuscule::logMemory([=](){ printf("+++ %p, size = %zu\n",ptr, count); });

  return ptr;
}
void* operator new  ( std::size_t count, std::align_val_t al) {

  // not available on osx:
  // auto ptr = aligned_alloc(count, al);

  void*ptr;
  if (int rc = posix_memalign(&ptr, static_cast<std::underlying_type_t<std::align_val_t>>(al), count)) {
      ptr = nullptr;
  }

  // make sure the capure is <= 16bytes:
  int c = count;
  int ali = static_cast<std::underlying_type_t<decltype(al)>>(al);
  imajuscule::logMemory([ptr, c, ali](){ printf("+++ %p, size = %d, alignment = %d\n",ptr, c, ali); });

  return ptr;
}
void operator delete(void* ptr) noexcept
{
  imajuscule::logMemory([=](){ printf("--- %p\n", ptr); });

  std::free(ptr);
}
void operator delete(void* ptr, std::align_val_t al) noexcept
{
  imajuscule::logMemory([=](){ printf("--- %p, alignment %lu\n", ptr, al); });

  std::free(ptr);
}
#  endif // IMJ_LOG_MEMORY
