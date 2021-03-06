

LIBRARY ieee;
USE ieee.std_logic_1164.ALL;
LIBRARY XilinxCoreLib;
ENTITY Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A IS
  PORT (
    clk : in std_logic := '0';
    a : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
    b : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
    ce : IN STD_LOGIC;
    sclr : IN STD_LOGIC;
    p : OUT STD_LOGIC_VECTOR(63 DOWNTO 0)
  );
END Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A;

ARCHITECTURE Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A_a OF Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A IS
COMPONENT wrapped_Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A
  PORT (
    clk : IN STD_LOGIC;
    a : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
    b : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
    ce : IN STD_LOGIC;
    sclr : IN STD_LOGIC;
    p : OUT STD_LOGIC_VECTOR(63 DOWNTO 0)
  );
END COMPONENT;

  FOR ALL : wrapped_Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A USE ENTITY XilinxCoreLib.mult_gen_v11_2(behavioral)
    GENERIC MAP (
      c_a_type => 0,
      c_a_width => 32,
      c_b_type => 0,
      c_b_value => "10000001",
      c_b_width => 32,
      c_ccm_imp => 0,
      c_ce_overrides_sclr => 0,
      c_has_ce => 1,
      c_has_sclr => 1,
      c_has_zero_detect => 0,
      c_latency => 8,
      c_model_type => 0,
      c_mult_type => 1,
      c_optimize_goal => 1,
      c_out_high => 63,
      c_out_low => 0,
      c_round_output => 0,
      c_round_pt => 0,
      c_verbosity => 0,
      c_xdevicefamily => "spartan6"
    );
BEGIN
U0 : wrapped_Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A
  PORT MAP (
    clk => clk,
    a => a,
    b => b,
    ce => ce,
    sclr => sclr,
    p => p
  );

END Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A_a;

configuration conf_87D7BDF9F83143F1B5E978F9A992495A of Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A is
  for Multiplier_1_87D7BDF9F83143F1B5E978F9A992495A_a end for; 
end conf_87D7BDF9F83143F1B5E978F9A992495A; 