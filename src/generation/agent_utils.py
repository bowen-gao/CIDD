"""
Agent utility classes for CIDD framework.

This module provides base classes for implementing LLM-based agents
that can collaborate in the molecular generation pipeline.
"""

import os
import requests
import json
from typing import List, Optional, Dict, Any


def send_request(context: List[Dict[str, str]], message: str) -> Optional[str]:
    """
    Send a request to OpenAI API.
    
    Args:
        context: List of message dictionaries for conversation context
        message: The message to send
        
    Returns:
        Response content from the API, or None if request failed
    """
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("OPENAI_API_KEY environment variable not set")
    
    # Use configurable API endpoint
    api_base = os.environ.get("OPENAI_API_BASE", "https://api.openai.com/v1")
    model = os.environ.get("OPENAI_MODEL", "gpt-4")
    
    url = f"{api_base}/chat/completions"
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    
    context.append({"role": "user", "content": message})
    payload = {
        "model": model,
        "messages": context,
        "temperature": 0.7
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers)
        response.raise_for_status()
        
        response_data = response.json()
        assistant_message = response_data['choices'][0]['message']['content']
        context.append({"role": "assistant", "content": assistant_message})
        return assistant_message
        
    except requests.exceptions.RequestException as e:
        print(f"API request failed: {e}")
        return None


class Agent:
    """
    Base class for LLM-based agents.
    
    Each agent maintains its own conversation context and can interact
    with language models to perform specific tasks.
    """
    
    def __init__(self, name: str):
        """
        Initialize the agent.
        
        Args:
            name: Name identifier for the agent
        """
        self.name = name
        self.context: List[Dict[str, str]] = []
    
    def send_request(self, context: List[Dict[str, str]], message: str) -> Optional[str]:
        """Send a request using the global send_request function."""
        return send_request(context, message)
    
    def execute(self, message: str) -> Optional[str]:
        """
        Execute a task by sending a message to the LLM.
        
        Args:
            message: The task description or question
            
        Returns:
            Response from the LLM, or None if failed
        """
        response = self.send_request(self.context, message)
        return response
    
    def reset_context(self):
        """Clear the conversation context."""
        self.context = []


class Part:
    """
    Container for multiple agents that work together.
    
    A Part represents a logical grouping of agents that collaborate
    on related tasks within the molecular generation pipeline.
    """
    
    def __init__(self, part_name: str):
        """
        Initialize the part.
        
        Args:
            part_name: Name identifier for the part
        """
        self.part_name = part_name
        self.agents: List[Agent] = []
    
    def add_agent(self, agent: Agent):
        """
        Add an agent to this part.
        
        Args:
            agent: The agent to add
        """
        self.agents.append(agent)
    
    def execute_all(self, message: str = "Executing coordinated task"):
        """
        Execute all agents in this part.
        
        Args:
            message: Task description to send to all agents
        """
        print(f"Executing Part: {self.part_name}")
        for agent in self.agents:
            agent.execute(message)
    
    def execute_agent(self, agent_name: str, message: str = "Executing specific task"):
        """
        Execute a specific agent by name.
        
        Args:
            agent_name: Name of the agent to execute
            message: Task description to send to the agent
        """
        for agent in self.agents:
            if agent.name == agent_name:
                agent.execute(message)
                break
        else:
            print(f"Agent '{agent_name}' not found in part '{self.part_name}'")
