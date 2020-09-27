
TEST(function, func) {
  folly::Function<void(void)> f;
  ASSERT_FALSE(static_cast<bool>(f));

  f = [](){
  };

  ASSERT_TRUE(static_cast<bool>(f));
}
