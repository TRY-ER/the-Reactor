import os
import pytest
import asyncio
import sys
import dotenv
dotenv.load_dotenv()

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from app.agents.one_shot_text_agent.agent import BaseOpenAIAgentInstance, DummyTextResponse, run_agent


class TestDummyTextResponseSchema:
    """Tests specifically for the DummyTextResponse output schema."""
    
    def test_dummy_text_response_schema_class_exists(self):
        """Test that DummyTextResponse class can be imported and is properly defined."""
        # Test that the class exists and has the expected attributes
        assert hasattr(DummyTextResponse, '__annotations__')
        assert 'response' in DummyTextResponse.__annotations__
        assert DummyTextResponse.__annotations__['response'] == str
    
    @pytest.mark.asyncio
    async def test_agent_with_schema_returns_structured_response(self):
        """Test that agent returns properly structured DummyTextResponse."""
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
        
        agent = BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq",
            output_schema=DummyTextResponse
        )
        
        response = await agent.run("Say hello in a friendly way")
        
        # Verify the response is of the expected type
        assert isinstance(response, DummyTextResponse)
        assert hasattr(response, 'response')
        assert isinstance(response.response, str)
        assert len(response.response) > 0
        
        # Check that the response contains meaningful content
        assert "hello" in response.response.lower() or "hi" in response.response.lower()
    
    def test_run_agent_function_compatibility(self):
        """Test that the run_agent function still works with the new schema."""
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
        
        # Test the convenience function
        response = run_agent("What is 1+1?")
        
        # The run_agent function should return the response
        # Since it uses the default agent which has DummyTextResponse schema
        assert response is not None
        
        # Check if it's a DummyTextResponse or string
        if isinstance(response, DummyTextResponse):
            assert hasattr(response, 'response')
            assert isinstance(response.response, str)
            assert len(response.response) > 0
        else:
            # Fallback check for string response
            assert isinstance(response, str)
            assert len(response) > 0
    
    @pytest.mark.asyncio 
    async def test_multiple_queries_consistent_format(self):
        """Test that multiple queries return consistent format."""
        if not os.getenv("GROQ_API_KEY"):
            pytest.skip("GROQ_API_KEY environment variable not set")
        
        agent = BaseOpenAIAgentInstance(
            model_name="groq/openai/gpt-oss-120b",
            endpoint_type="groq",
            output_schema=DummyTextResponse
        )
        
        queries = [
            "What is the capital of France?",
            "Tell me a short joke",
            "Explain what AI is in one sentence"
        ]
        
        for query in queries:
            response = await agent.run(query)
            
            # Each response should be consistently formatted
            assert isinstance(response, DummyTextResponse)
            assert hasattr(response, 'response')
            assert isinstance(response.response, str)
            assert len(response.response.strip()) > 0
            
            # Basic content validation
            assert any(c.isalnum() for c in response.response)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])