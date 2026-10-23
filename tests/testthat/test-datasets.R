test_that("all practice datasets load with expected columns and non-empty rows", {
  # 1. wage_data
  data(wage_data, package = "phoring")
  expect_equal(nrow(wage_data), 600)
  expect_true(all(c("wage", "educ", "exper", "tenure") %in% names(wage_data)))
  
  # 2. house_prices
  data(house_prices, package = "phoring")
  expect_equal(nrow(house_prices), 506)
  expect_true(all(c("price", "rooms", "crime", "nox") %in% names(house_prices)))
  
  # 3. college_distance
  data(college_distance, package = "phoring")
  expect_equal(nrow(college_distance), 800)
  expect_true(all(c("wage", "educ", "near_college") %in% names(college_distance)))
  
  # 4. firm_investment
  data(firm_investment, package = "phoring")
  expect_equal(nrow(firm_investment), 200)
  expect_true(all(c("firm_id", "year", "inv", "value", "capital") %in% names(firm_investment)))
  
  # 5. crime_panel
  data(crime_panel, package = "phoring")
  expect_equal(nrow(crime_panel), 200)
  expect_true(all(c("county_id", "year", "crime_rate", "police") %in% names(crime_panel)))
  
  # 6. labor_supply
  data(labor_supply, package = "phoring")
  expect_equal(nrow(labor_supply), 753)
  expect_true(all(c("inlf", "nwifeinc", "educ", "kidslt6") %in% names(labor_supply)))
  
  # 7. min_wage
  data(min_wage, package = "phoring")
  expect_equal(nrow(min_wage), 400)
  expect_true(all(c("store_id", "treated", "post", "fte") %in% names(min_wage)))
  
  # 8. macro_data
  data(macro_data, package = "phoring")
  expect_equal(nrow(macro_data), 120)
  expect_true(all(c("quarter", "inflation", "unemployment") %in% names(macro_data)))
})
