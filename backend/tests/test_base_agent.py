import os
import pytest
import asyncio
import sys
import dotenv
dotenv.load_dotenv()

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from app.agents.one_shot_text_agent.agent import BaseOpenAIAgentInstance, DummyTextResponse


class TestBaseOpenAIAgentInstance:
    """Simple tests for BaseOpenAIAgentInstance to check if methods return valid string outputs."""
    
    @pytest.fixture
    def agent_instance(self):
        """Create a BaseOpenAIAgentInstance for testing."""
        # Skip test if no API key is available
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
        
        return BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq",
            agent_name="test_agent",
            instructions="Test agent for unit testing purposes."
        )
    
    @pytest.mark.asyncio
    async def test_run_returns_valid_string(self, agent_instance):
        """Test that run method returns a valid string response."""
        response = await agent_instance.run("Hello")
        
        # Check if response is a string and not empty
        assert isinstance(response, str)
        assert len(response) > 0
        assert response.strip() != ""
    
    @pytest.mark.asyncio
    async def test_run_returns_dummy_text_response_format(self, agent_instance):
        """Test that run method returns response in DummyTextResponse format when using output schema."""
        # Test with explicit output schema
        agent_with_schema = BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq",
            agent_name="test_agent_schema",
            instructions="Test agent that should return structured response.",
            output_schema=DummyTextResponse
        )
        
        response = await agent_with_schema.run("Say hello")
        
        # The response should be a DummyTextResponse instance
        assert isinstance(response, DummyTextResponse)
        assert hasattr(response, 'response')
        assert isinstance(response.response, str)
        assert len(response.response) > 0
        assert response.response.strip() != ""
    
    @pytest.mark.asyncio
    async def test_run_response_content_validation(self, agent_instance):
        """Test that the response content is valid and meaningful."""
        response = await agent_instance.run("What is 2+2?")
        
        # Extract text from response regardless of format
        if isinstance(response, DummyTextResponse):
            response_text = response.response
        else:
            response_text = str(response)
        
        # Check if response is a string and contains meaningful content
        assert isinstance(response_text, str)
        assert len(response_text) > 0
        assert response_text.strip() != ""
        # Basic check that it's not just whitespace or empty
        assert any(c.isalnum() for c in response_text)
    
    @pytest.mark.asyncio
    async def test_run_streamed_returns_valid_responses(self, agent_instance):
        """Test that run_streamed method returns valid responses."""
        responses = []
        
        async for response in agent_instance.run_streamed("Say hello"):
            # Just collect a few responses to test
            responses.append(response)
            if len(responses) >= 3:  # Just test first few responses
                break
        
        # Check that we got some responses
        assert len(responses) > 0
        
        # Each response should be a valid stream event
        for response in responses:
            assert response is not None
    
    @pytest.mark.asyncio
    async def test_output_schema_consistency(self, agent_instance):
        """Test that the agent consistently uses the output schema."""
        # Test multiple queries to ensure consistent output format
        queries = ["Hello", "What is AI?", "Tell me a joke"]
        
        for query in queries:
            response = await agent_instance.run(query)
            
            # Check the type based on the agent's output_schema
            if agent_instance.output_schema == DummyTextResponse:
                assert isinstance(response, DummyTextResponse)
                assert hasattr(response, 'response')
                assert isinstance(response.response, str)
                assert len(response.response) > 0
            else:
                # When no schema is set, agent returns plain string
                assert isinstance(response, str)
                assert len(response) > 0
    
    def test_get_model_returns_something(self, agent_instance):
        """Test that get_model method returns something (not None)."""
        model = agent_instance.get_model()
        
        # Just check that it returns something
        assert model is not None
    
    def test_output_schema_initialization(self):
        """Test that the agent is properly initialized with the output schema."""
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
            
        agent = BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq",
            output_schema=DummyTextResponse
        )
        
        # Check that the output schema is properly set
        assert agent.output_schema == DummyTextResponse
        assert agent.agent.output_type == DummyTextResponse
    
    def test_default_output_schema(self):
        """Test that the agent has no default output schema."""
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
            
        agent = BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq"
        )
        
        # Check that the default output schema is None (no structured output by default)
        assert agent.output_schema is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
