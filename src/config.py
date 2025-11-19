"""
CIDD Configuration Module

This module contains configuration settings and environment variables
for the CIDD framework. Modify these settings to match your environment.
"""

import os
from typing import Optional


class CIDDConfig:
    """Configuration class for CIDD framework."""
    
    # API Configuration
    OPENAI_API_KEY: Optional[str] = os.getenv('OPENAI_API_KEY')
    OPENAI_API_BASE: str = os.getenv('OPENAI_API_BASE', 'https://api.openai.com/v1')
    OPENAI_MODEL: str = os.getenv('OPENAI_MODEL', 'gpt-4')
    
    # Volcengine API (optional)
    VOLCENGINE_API_KEY: Optional[str] = os.getenv('VOLCENGINE_API_KEY')
    
    # Data Paths
    DATA_ROOT: str = os.getenv('CIDD_DATA_ROOT', './data')
    OUTPUT_DIR: str = os.getenv('CIDD_OUTPUT_DIR', './results')
    
    # Docking Configuration
    VINA_EXHAUSTIVENESS: int = int(os.getenv('VINA_EXHAUSTIVENESS', '16'))
    DOCKING_TIMEOUT: int = int(os.getenv('DOCKING_TIMEOUT', '300'))
    
    # Generation Parameters
    MAX_ITERATIONS: int = int(os.getenv('CIDD_MAX_ITERATIONS', '5'))
    SCAFFOLD_EXPLORE: int = int(os.getenv('CIDD_SCAFFOLD_EXPLORE', '10'))
    SIDECHAIN_EXPLORE: int = int(os.getenv('CIDD_SIDECHAIN_EXPLORE', '3'))
    
    # Parallel Processing
    MAX_WORKERS: int = int(os.getenv('CIDD_MAX_WORKERS', '10'))
    TASK_TIMEOUT: int = int(os.getenv('CIDD_TASK_TIMEOUT', '10800'))
    
    @classmethod
    def validate_config(cls) -> bool:
        """
        Validate configuration settings.
        
        Returns:
            True if configuration is valid, False otherwise
        """
        issues = []
        
        if not cls.OPENAI_API_KEY:
            issues.append("OPENAI_API_KEY not set - LLM features will not work")
        
        if not os.path.exists(cls.DATA_ROOT):
            issues.append(f"Data root directory does not exist: {cls.DATA_ROOT}")
        
        if issues:
            print("⚠️  Configuration Issues:")
            for issue in issues:
                print(f"   - {issue}")
            return False
        
        return True
    
    @classmethod
    def print_config(cls):
        """Print current configuration settings."""
        print("🔧 CIDD Configuration:")
        print(f"   Data Root: {cls.DATA_ROOT}")
        print(f"   Output Dir: {cls.OUTPUT_DIR}")
        print(f"   OpenAI Model: {cls.OPENAI_MODEL}")
        print(f"   Max Workers: {cls.MAX_WORKERS}")
        print(f"   Vina Exhaustiveness: {cls.VINA_EXHAUSTIVENESS}")
        print(f"   API Key Set: {'✓' if cls.OPENAI_API_KEY else '✗'}")


# Global configuration instance
config = CIDDConfig()