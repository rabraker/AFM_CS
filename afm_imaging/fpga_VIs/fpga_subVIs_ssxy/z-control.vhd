-- Description : A RAM based FIFO buffer
--
--

library ieee;
-- use ieee.fixed_pkg.all;

use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity z_control is
  port (
    clk         : in  std_logic;                            -- clock
    sclr        : in  std_logic;
    ce          : in  std_logic;
    state       : in std_logic_vector(7 downto 0);
    u_kmin1     : in std_logic_vector(15 downto 0);
    u_k         : out std_logic_vector(15 downto 0);
    fast_thresh : in std_logic_vector(15 downto 0);
    z_err       : in std_logic_vector(15 downto 0);
    K_I         : in std_logic_vector(15 downto 0);
    r_fast      : in std_logic_vector(15 downto 0);
    r_slow      : in std_logic_vector(15 downto 0)
    );

end entity z_control;


architecture synth of z_control is

  -- signal fifo_state_cur_reg  : fifo_state := fifo_state_init;  -- current fifo state
  -- signal fifo_state_next_reg : fifo_state := fifo_state_init;  -- next fifo state
  -- signal addr_diff_reg       : std_logic_vector(AWIDTH downto 0) := (others => '0');
  -- signal RAM_we              : slbit := '0';

begin  -- architecture synth

  -- proc_clk: process(all)
  -- begin
  --   -- if rising_edge(clk) then
  --   --   fifo_state_cur_reg <= fifo_state_next_reg;
  --   -- end if;
  -- end process proc_clk;


  proc_zcontrol: process(clk)
  begin
    if rising_edge(clk) then
      if ce = '1' then
      if sclr = '1' then
        u_k <= (others => '0');
      elsif (unsigned(state) = 1 or unsigned(state) = 2) and u_kmin1 > fast_thresh then
        u_k <= std_logic_vector(signed(u_kmin1) + signed(r_fast));
      -- elsif unsigned(state) = 1  then
      --   u_k <= u_kmin1;
      -- elsif unsigned(state) = 2 then
      --   u_k <= std_logic_vector(signed(u_kmin1) - signed(r_slow));
      -- elsif unsigned(state) = 3 then
      --   u_k <= std_logic_vector(signed(u_kmin1) + signed(K_I)*signed(z_err))
      end if;
      end if;
  end if;


  end process proc_zcontrol;

end architecture synth;
